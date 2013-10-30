package CracTools::Annotator;

use strict;
use warnings;

use Carp;
use Data::Dumper;
use CracTools::GFF::Annotation;
use CracTools::GFF::Query2;
use CracTools::Const;

our @type_prot = ('3PRIM_UTR_sense','CDS_sense','5PRIM_UTR_sense','INXON_sense','INTRON_sense','3PRIM_UTR_antisense','CDS_antisense','5PRIM_UTR_antisense','INXON_antisense','INTRON_antisense','INTER_PROXIMAL','INTER_DISTAL_EST','INTER_DISTAL');

our @type_noncoding = ("small_ncRNA","lincRNA","other_lncRNA","other_noncodingRNA");

our %annotation_priority = ($type_prot[0] => 1,
                            $type_prot[1] => 2,
                            $type_prot[2] => 3,
                            $type_prot[3] => 4,
                            $type_prot[4] => 5,
                            $type_prot[5] => 6,
                            $type_prot[6] => 7,
                            $type_prot[7] => 8,
                            $type_prot[8] => 9,
                            $type_prot[9] => 10,
                            $type_prot[10] => 11,
                            $type_prot[11] => 12,
                            $type_prot[12] => 13);

our %priority_annotation = reverse %annotation_priority;


=head1 METHODS

=head2 new

  Arg [1] : String - $gff_file
            GFF file to perform annotation

  Example     : my $annotation = CracTools::GFF::Annotation->new($gff_line);
  Description : Create a new CracTools::GFF::Annotation object
                If a gff line is passed in argument, the line will be parsed
                and loaded.
  ReturnType  : CracTools::GFF::Query
  Exceptions  : none

=cut

sub new {
  my $class = shift;
  my $gff_file = shift;

  if(!defined $gff_file) {
    croak "Missing GFF file argument in CracTools::Annotator constructor";
  }

  my $self = bless {
    gff_file => $gff_file,
  }, $class;

  $self->_init();

  return $self;
}

sub foundSameGene {
  my $self = shift;
  my ($chr,$pos_start1,$pos_end1,$pos_start2,$pos_end2,$strand) = @_;
  my ($candidates1_ref,$genes1_ref) = $self->_getAnnotationCandidates($chr,$pos_start1,$pos_end1,$strand);
  my ($candidates2_ref,$genes2_ref) = $self->_getAnnotationCandidates($chr,$pos_start2,$pos_end2,$strand);
  my $found_same_gene = 0;
  foreach my $gene_key (keys %{$genes1_ref}) {
    if(defined $genes2_ref->{$gene_key}) {
      $found_same_gene = 1;
      last;
    }
  }
  return $found_same_gene;
}

# TODO add error handles in case GFF3 doesnt have required informations
sub getAnnotation {
  my $self = shift;
  my ($chr,$pos_start,$pos_end,$strand) = @_;

  my %annotation;

  my ($candidates_ref,$genes_ref) = $self->_getAnnotationCandidates($chr,$pos_start,$pos_end,$strand);

  # First we take a look to annotations on the original strand
  my ($best_coding_candidate,$best_coding_priority,$non_coding) = (undef,13,undef);
  foreach my $candidate (values %{$candidates_ref}) {
    my $flags = $candidate->{flags};
    if(!defined $candidate->{mRNA}) {
      print Dumper($candidate);
      #print $candidate->{exon}->attribute('ID'),"\n";
    }
    if($candidate->{mRNA}->attribute('type') =~ /protein_coding/i) {
      my $candidate_type = $self->_getCandidateType($candidate).'_sense';
      if($annotation_priority{$candidate_type} < $best_coding_priority) {
        $best_coding_candidate = $candidate;
        $best_coding_priority = $annotation_priority{$candidate_type};
      }
    } elsif(!defined $non_coding) {
      $non_coding = $candidate;
    }
  }

  # If we haven't found a coding candidate on the original strand
  # we look to the opposite strand (only for coding candidates)
  if(!defined $best_coding_candidate) {
    my ($candidates_opposite_strand_ref,$genes_opposite_strand_ref) = $self->_getAnnotationCandidates($chr,$pos_start,$pos_end,$strand*-1);

    # Merge genes hashes
    # We have to keep data of the previous one if we have found
    # A non-coding.
    @{$genes_ref}{keys %{$genes_opposite_strand_ref}} = values %{$genes_opposite_strand_ref};

    foreach my $candidate (values %{$candidates_opposite_strand_ref}) {
      my $flags = $candidate->{flags};
      if($candidate->{mRNA}->attribute('type') =~ /protein_coding/i) {
        my $candidate_type = $self->_getCandidateType($candidate).'_antisense';
        if($annotation_priority{$candidate_type} < $best_coding_priority) {
          $best_coding_candidate = $candidate;
          $best_coding_priority = $annotation_priority{$candidate_type};
        }
      }
    }
  }

  # If we haven't found a coding candidate neither on the original
  # strand nor on the opposite strand, we try looking with a wider
  # window.
  if(!defined $best_coding_candidate && $self->neighborhoodSearch) {
    my $distance_max = $self->distanceMax();
    my $intergenic_threshold = $self->intergenicThreshold();

    # First we look on the 5PRIM_UTR side
    my ($candidates_5prim_utr_ref,$genes_5prim_utr_ref) = $self->_getAnnotationCandidates($chr,$pos_start - $distance_max,$pos_start-1,$strand);

    my $best_5prim_utr_candidate_distance = $distance_max;
    my $best_5prim_utr_candidate;

    foreach my $candidate (values %{$candidates_5prim_utr_ref}) {
      my $mRNA = $candidate->{mRNA};
      # We have found a closer gene
      my $gene_distance = $pos_start - $mRNA->end;
      if($gene_distance < $best_5prim_utr_candidate_distance && $gene_distance > 0) {
        $best_5prim_utr_candidate = $candidate;
        $best_5prim_utr_candidate_distance = $gene_distance; 
      }
    }

    # Then we look on the 3PRIM_UTR side
    my ($candidates_3prim_utr_ref,$genes_3prim_utr_ref) = $self->_getAnnotationCandidates($chr,$pos_end+1,$pos_end+$distance_max,$strand);

    my $best_3prim_utr_candidate_distance = $distance_max;
    my $best_3prim_utr_candidate;

    foreach my $candidate (values %{$candidates_3prim_utr_ref}) {
      my $mRNA = $candidate->{mRNA};
      # We have found a closer gene
      my $gene_distance = $mRNA->start - $pos_end;
      if($gene_distance < $best_3prim_utr_candidate_distance && $gene_distance > 0) {
        $best_3prim_utr_candidate = $candidate;
        $best_3prim_utr_candidate_distance = $gene_distance; 
      }
    }

    # If we are on the reverse strand we swap 5prim_utr and 3prim_utr
    # candidates.
    if($strand == -1) {
      ($best_3prim_utr_candidate_distance,$best_5prim_utr_candidate_distance) = ($best_5prim_utr_candidate_distance,$best_3prim_utr_candidate_distance);
      ($best_3prim_utr_candidate,$best_5prim_utr_candidate) = ($best_5prim_utr_candidate,$best_3prim_utr_candidate);
      ($genes_3prim_utr_ref,$genes_5prim_utr_ref) = ($genes_5prim_utr_ref,$genes_3prim_utr_ref);
    }

    # Now we fill the annotation hash with the bests 3prim_utr
    # and 5prim_utr candidates
    if(defined $best_3prim_utr_candidate) {
      my $mRNA = $best_3prim_utr_candidate->{mRNA};
      my $gene = $genes_3prim_utr_ref->{$mRNA->attribute('Parent')};
      $annotation{hugo_3prim} = $gene->attribute('Name');
      $annotation{desc_3prim} = $mRNA->attribute('type');
      $annotation{id_3prim} = $gene->attribute('ID');
      $annotation{distance_of_3prim_gene} = $best_3prim_utr_candidate_distance;
    }

    if(defined $best_5prim_utr_candidate) {
      my $mRNA = $best_5prim_utr_candidate->{mRNA};
      my $gene = $genes_5prim_utr_ref->{$mRNA->attribute('Parent')};
      $annotation{hugo_5prim} = $gene->attribute('Name');
      $annotation{desc_5prim} = $mRNA->attribute('type');
      $annotation{id_5prim} = $gene->attribute('ID');
      $annotation{distance_of_5prim_gene} = $best_5prim_utr_candidate_distance;
    }

    # If we are close (according to the 'intergenic_threshold) to the 3prim UTR
    # of a gene we have an INTER_PROXIMAL priority.
    # Either we are in an intergenic zone.
    if(defined $best_5prim_utr_candidate && $best_5prim_utr_candidate_distance <= $intergenic_threshold) {
      $best_coding_priority = $annotation_priority{INTER_PROXIMAL};
    } else {
      $best_coding_priority = $annotation_priority{INTER_DISTAL};
    }
  }

  # We fill the annotation hash with best candidates data
  $annotation{annot} = $priority_annotation{$best_coding_priority};
  $annotation{priority} = $best_coding_priority;
  if(defined $best_coding_candidate) {
    my $mRNA = $best_coding_candidate->{mRNA};
    my $gene = $genes_ref->{$mRNA->attribute('Parent')};
    my $exon =  $best_coding_candidate->{exon};
    $annotation{hugo} = $gene->attribute('Name');
    $annotation{id} = $gene->attribute('ID');
    $annotation{desc} = $mRNA->attribute('type');
    if(defined $exon) {
      $annotation{exon} = $exon;
    }
  }
  if(defined $non_coding) {
    my $mRNA_non_coding = $non_coding->{mRNA};
    my $gene_non_coding = $genes_ref->{$mRNA_non_coding->attribute('Parent')};
    $annotation{hugo_non_coding} = $gene_non_coding->attribute('Name');
    $annotation{id_non_coding} = $gene_non_coding->attribute('ID');
    $annotation{desc_non_coding} = $mRNA_non_coding->attribute('type');
  }

  # TEST printing annotation hash key/values
  #while(my ($k,$v) = each %annotation) {
  #  print "$k => $v\n";
  #}

  # Return a reference to the annotation hash
  return \%annotation;
}

sub distanceMax {
  my $self = shift;
  my $distance_max = shift;
  if(defined $distance_max) {
    $self->{DISTANCE_MAX} = $distance_max;
  }elsif(!defined $self->{DISTANCE_MAX}) {
    $self->{DISTANCE_MAX} = $CracTools::Const::ANNOTATION_DISTANCE_MAX;
  }
  return $self->{DISTANCE_MAX};
}  

sub intergenicThreshold {
  my $self = shift;
  my $intergenic_threshold = shift;
  if(defined $intergenic_threshold) {
    $self->{INTERGENIC_THRESHOLD} = $intergenic_threshold;
  }elsif(!defined $self->{INTERGENIC_THRESHOLD}) {
    $self->{INTERGENIC_THRESHOLD} = $CracTools::Const::INTERGENIC_THRESHOLD;
  }
  return $self->{INTERGENIC_THRESHOLD};
}  

sub neighborhoodSearch {
  my $self = shift;
  my $neighborhood_search = shift;
  if(defined $neighborhood_search) {
    $self->{NEIGHBORHOOD_SEARCH} = $neighborhood_search;
  }elsif(!defined $self->{NEIGHBORHOOD_SEARCH}) {
    $self->{NEIGHBORHOOD_SEARCH} = $CracTools::Const::NEIGHBORHOOD_SEARCH;
  }
  return $self->{NEIGHBORHOOD_SEARCH};
}

sub _init {
  my $self = shift;

  # Create a GFF file to query exons
  my $gff_query = CracTools::GFF::Query2->new($self->{gff_file});
  $self->{gff_query} = $gff_query;

}


sub _getAnnotationCandidates {
  my $self = shift;
  my ($chr,$pos_start,$pos_end,$strand) = @_;

  # get GFF annotations that overlap the region to annotate
  my $annotations = $self->{gff_query}->fetchByRegion($chr,$pos_start,$pos_end,$strand);

  my %candidates;
  my %genes;

  foreach my $annot_line (@{$annotations}) {
    my $annot = CracTools::GFF::Annotation->new($annot_line,'gff3');
    my $annot_priority = 0; 
    if(defined $annot->feature) {
      if ($annot->feature =~ /exon/i) {
        # An exon can have multiple Parents
        my @parents = $annot->parents();
        foreach my $parent (@parents) {
          $candidates{$parent}{exon} = $annot;
          if($annot->start <= $pos_start && $annot->end >= $pos_end) {
            $candidates{$parent}{flags} += 1;
          } else {
            $candidates{$parent}{flags} += 2;
          }
        }
      } elsif($annot->feature =~ /five/i) {
        $candidates{$annot->attribute('Parent')}{flags} += 4; 
        $candidates{$annot->attribute('Parent')}{five} = $annot;
      } elsif($annot->feature =~ /three/i) {
        $candidates{$annot->attribute('Parent')}{three} = $annot;
        $candidates{$annot->attribute('Parent')}{flags} += 8; 
      } elsif($annot->feature =~ /cds/i) {
        $candidates{$annot->attribute('Parent')}{cds} = $annot;
        $candidates{$annot->attribute('Parent')}{flags} += 16; 
      } elsif($annot->feature =~ /mRNA/i) {
        $candidates{$annot->attribute('ID')}{flags} += 32; 
        $candidates{$annot->attribute('ID')}{mRNA} = $annot; 
      } elsif($annot->feature =~ /gene/i) {
        $genes{$annot->attribute('ID')} = $annot; 
      }
    }
  }
  return (\%candidates,\%genes);
}

sub _getCandidateType {
  my $self = shift;
  my $candidate_ref = shift;
  my $flags = $candidate_ref->{flags};
  my $candidate_type;
  if ($flags & 1) {
    if ($flags & 16) {
      $candidate_type = 'CDS';
    } elsif ($flags & 4) {
      $candidate_type = '5PRIM_UTR';
    } elsif ($flags & 8) {
      $candidate_type = '3PRIM_UTR';
    } else {
      # TODO what do we do?
      # CDS is just for now
      $candidate_type = 'CDS';
    }
  } elsif ($flags & 2) {
    $candidate_type = 'INXON';
  } else {
    $candidate_type = 'INTRON';
  }
  return $candidate_type;
}
