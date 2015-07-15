package CracTools::Annotator;
# ABSTRACT: Generic annotation base on CracTools::GFF::Query::File

use strict;
use warnings;

use Carp;
use CracTools::GFF::Annotation;
use CracTools::Interval::Query;
use CracTools::Interval::Query::File;
use List::Util qw[min max];
use CracTools::Const;

=head1 SYNOPSYS

0-based coordinate system
closed [a,b] intervals

=head1 METHODS

=head2 new

  Arg [1] : String - $gff_file
            GFF file to perform annotation
  Arg [2] : String - $mode
            "fast", "light"

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
  my $mode = shift;

  $mode = 'light' if !defined $mode;

  if(!defined $gff_file) {
    croak "Missing GFF file argument in CracTools::Annotator constructor";
  }

  my $self = bless {
    gff_file => $gff_file,
    mode => $mode,
  }, $class;

  $self->_init();

  return $self;
}

=head2 mode 

=cut 

sub mode {
  my $self = shift;
  return $self->{mode};
}

=head2 foundGene

  Arg [1] : String - chr
  Arg [2] : String - pos_start
  Arg [3] : String - pos_end
  Arg [4] : String - strand

  Description : Return true if there is an exon of a gene is this interval
  ReturnType  : Boolean
  Exceptions  : none

=cut

sub foundGene {
  my $self = shift;
  my ($chr,$pos_start,$pos_end,$strand) = @_;
  my @candidates = @{ $self->getAnnotationCandidates($chr,$pos_start,$pos_end,$strand)};
  return (scalar @candidates > 0);
}

=head2 foundSameGene

  Arg [1] : String - chr
  Arg [2] : String - pos_start1
  Arg [3] : String - pos_end1
  Arg [4] : String - pos_start2
  Arg [5] : String - pos_end1
  Arg [6] : String - strand

  Description : Return true if a gene is the same gene is found is the two intervals.
  ReturnType  : Boolean
  Exceptions  : none

=cut

sub foundSameGene {
  my $self = shift;
  my ($chr,$pos_start1,$pos_end1,$pos_start2,$pos_end2,$strand) = @_;
  my @candidates1 = @{ $self->getAnnotationCandidates($chr,$pos_start1,$pos_end1,$strand)};
  my @candidates2 = @{ $self->getAnnotationCandidates($chr,$pos_start2,$pos_end2,$strand)};
  my $found_same_gene = 0;
  my @genes1;
  my @genes2;
  foreach my $candi1 (@candidates1) {
    if(defined $candi1->{gene}) {
      push @genes1,$candi1->{gene}->attribute('ID');
    }
  }
  foreach my $candi2 (@candidates2) {
    if(defined $candi2->{gene}) {
      push @genes2,$candi2->{gene}->attribute('ID');
    }
  }
  foreach my $gene_id (@genes1) {
    foreach (@genes2) {
      if($gene_id eq $_) {
        $found_same_gene = 1;
        last;
      }
    }
    last if $found_same_gene == 1;
  }
  return $found_same_gene;
}

=head2 getBestAnnotationCandidate

  Arg [1] : String - chr
  Arg [2] : String - pos_start
  Arg [3] : String - pos_end
  Arg [4] : String - strand
  Arg [5] : (Optional) Subroutine - see C<getCandidatePriorityDefault> for more details
  Arg [6] : (Optional) Subroutine - see C<compareTwoCandidatesDefault> for more details

  Description : Return best annotation candidate according to the priorities given
                by the subroutine in argument.
  ReturnType  : Hash( feature_name => CracTools::GFF::Annotation, ...), Int(priority), String(type)

=cut

sub getBestAnnotationCandidate {
  my $self = shift;
  my ($best_candidates,$best_priority,$best_type) = $self->getBestAnnotationCandidates(@_);
  if(@{$best_candidates}) {
    return $best_candidates->[0],$best_priority,$best_type;
  } else {
    return undef,undef,undef;
  }
}

=head2 getBestAnnotationCandidates

  Arg [1] : String - chr
  Arg [2] : String - pos_start
  Arg [3] : String - pos_end
  Arg [4] : String - strand
  Arg [5] : (Optional) Subroutine - see C<getCandidatePriorityDefault> for more details
  Arg [6] : (Optional) Subroutine - see C<compareTwoCandidatesDefault> for more details

  Description : Return best annotation candidates according to the priorities given
                by the subroutine(s) in argument.
  ReturnType  : ArrayRef[HashRef{ feature_name => CracTools::GFF::Annotation, ...}, ...], Int(priority), String(type)

=cut

sub getBestAnnotationCandidates {
  my $self = shift;
  my ($chr,$pos_start,$pos_end,$strand,$prioritySub,$compareSub) = @_;

  if(!defined $prioritySub && !defined $compareSub) {
    $prioritySub = \&getCandidatePriorityDefault unless defined $prioritySub;
    $compareSub = \&compareTwoCandidatesDefault unless defined $compareSub;
  }

  my @candidates = @{ $self->getAnnotationCandidates($chr,$pos_start,$pos_end,$strand)};
  my @best_candidates;
  my ($best_priority,$best_type);
  foreach my $candi (@candidates) {
    my ($priority,$type);
    ($priority,$type) = $prioritySub->($pos_start,$pos_end,$candi) if defined $prioritySub;
    if(defined $priority && $priority != -1) {
      if(!defined $best_priority) {
        $best_priority = $priority;
        push @best_candidates, $candi;
        $best_type = $type;
      } elsif($priority < $best_priority) {
        @best_candidates = ($candi);
        $best_priority = $priority;
        $best_type = $type;
      }
      #we should compare two candidates with equal priority to always choose the one
      elsif (!defined $priority || $priority == $best_priority){
        my $candidate_chosen;
        my $found_better_candidate = 0;
        foreach my $best_candidate (@best_candidates) {
          $candidate_chosen = $compareSub->($best_candidate,$candi,$pos_start,$pos_end) if defined $compareSub;
          # They are both equal
          if (!defined $candidate_chosen) {
            # We cannnot say if this candidate is better
            next;
          } elsif ($candidate_chosen == $candi) {
            # We have found a better candidate that previously register ones
            # we save it and remove the others
            @best_candidates = ($candi);
            $found_better_candidate = 1;
            last;
          } else {
            # The better candidate is not "candi", so this candidates
            # does not belong the the best_candidate array.
            # We can stop looping
            $found_better_candidate = 1;
            last;
          }
        }
        push @best_candidates, $candi if !$found_better_candidate;
      }
    }
  }
  # TODO We should not return variable in that order,
  # it is not easy to only retrieve the best candidatse...
  return \@best_candidates,$best_priority,$best_type;
}

=head2 getAnnotationCandidates

  Arg [1] : String - chr
  Arg [2] : String - pos_start
  Arg [3] : String - pos_end
  Arg [4] : String - strand

  Description : Return an array with all annotation candidates overlapping the
                chromosomic region.
  ReturnType  : Array of Hash( feature_name => CracTools::GFF::Annotation, ...)

=cut

sub getAnnotationCandidates {
  my $self = shift;
  my ($chr,$pos_start,$pos_end,$strand) = @_;
  # TODO if no strand is provided we should return annotations from both strands

  # get GFF annotations that overlap the region to annotate
  my $annotations = $self->{gff_query}->fetchByRegion($chr,$pos_start,$pos_end,$strand);
  # get a ref of an array of hash of candidates
  my $candidatates = $self->_constructCandidatesFromAnnotation($annotations);
  return $candidatates;
}

=head2 getAnnotationNearestDownCandidates

  Arg [1] : String - chr
  Arg [2] : String - pos_start
  Arg [3] : String - strand

  Description : Return an array with all annotation candidates nearest down the
                query region (without overlap).
  ReturnType  : Array of Hash( feature_name => CracTools::GFF::Annotation, ...)

=cut

sub getAnnotationNearestDownCandidates {
  my $self = shift;
  my ($chr,$pos_start,$strand) = @_;

  # get GFF annotations that overlap the pos_start to annotate
  my $annotations_overlap = $self->{gff_query}->fetchByLocation($chr,$pos_start,$strand);
  # get GFF annotations of nearest down intervals that not overlaped [pos_start,pos_end] pos 
  my @annotations_down;

  push @annotations_down, @{$self->{gff_query}->fetchAllNearestDown($chr,$pos_start,$strand)};

  # get a ref of an array of hash of candidates
  my @annotations = (@$annotations_overlap,@annotations_down);
  my $candidatates = $self->_constructCandidatesFromAnnotation(\@annotations);
  return $candidatates;
}

=head2 getAnnotationNearestUpCandidates

  Arg [1] : String - chr
  Arg [2] : String - pos_end
  Arg [3] : String - strand

  Description : Return an array with all annotation candidates nearest up the
                query region (without overlap).
  ReturnType  : ArrayRef of HashRef{ feature_name => CracTools::GFF::Annotation, ...}

=cut

sub getAnnotationNearestUpCandidates {
  my $self = shift;
  my ($chr,$pos_end,$strand) = @_;

  # get GFF annotations that overlap the pos_end to annotate
  my $annotations_overlap = $self->{gff_query}->fetchByLocation($chr,$pos_end,$strand);
  # get GFF annotations of nearest up intervals that not overlaped [pos_start,pos_end] pos 
  my @annotations_up;

  push @annotations_up, @{$self->{gff_query}->fetchAllNearestUp($chr,$pos_end,$strand)};

  # get a ref of an array of hash of candidates
  my @annotations = (@$annotations_overlap,@annotations_up);
  my $candidatates = $self->_constructCandidatesFromAnnotation(\@annotations);
  return $candidatates;
}

=head2 getCandidatePriorityDefault

  Arg [1] : String - pos_start
  Arg [2] : String - pos_end
  Arg [3] : hash - candidate

  Description : Default method used to give a priority to a candidate.
                You can create your own priority method to fit your specific need
                for selecting the best annotation.
                The best priority is 0. A priority of -1 means that this candidate
                should be avoided.
  ReturnType  : Array ($priority,$type) where $priority is an integer and $type a string

=cut

sub getCandidatePriorityDefault {
  my ($pos_start,$pos_end,$candidate) = @_;
  my ($priority,$type) = (-1,'');
  my ($mRNA,$exon) = ($candidate->{mRNA},$candidate->{exon});
  if(defined $mRNA) {
    if($mRNA->attribute('type') =~ /protein_coding/i) {
      if(defined $exon) {
        if(($exon->start <= $pos_start) && ($exon->end >= $pos_end)) {
          $priority = 1;
          if(defined $candidate->{three}) {
            $type = '3PRIM_UTR';
          } elsif(defined $candidate->{five}) {
            $type = '5PRIM_UTR';
          # } elsif(defined $candidate->{cds}) {
          #   $type = 'CDS';
          } else {
            $type = 'EXON';
          }
        } else {
          $priority = 2;
          $type = 'INXON';
        }
      } else {
	$priority = 4;
	$type = 'INTRON';
      }
    } else {
      if(defined $exon) {
        if(($exon->start <= $pos_start) && ($exon->end >= $pos_end)) {
          $priority = 3;
          $type = 'NON_CODING';
        }
      }
    }
  }
  return ($priority,$type);
}

=head2 compareTwoCandidatesDefault

  Arg [1] : hash - candidate1
  Arg [2] : hash - candidate2
  Arg [3] : pos_start (position start that has been queried)
  Arg [4] : pos_end (position end that has been queried)

  Description : Default method used to chose the best candidat when priority are equals
                You can create your own priority method to fit your specific need
                for selecting the best candidat.
  ReturnType  : hash - best candidat or undef if we cannot decide which candidate is the best

=cut
sub compareTwoCandidatesDefault{
  my ($candidate1,$candidate2,$pos_start) = @_;
  # If both candidates are exons we try to find wich one is closer to the pos_start of the region to annotate
  if ($candidate1->{exon} && $candidate2->{exon}) { 
    my $dist1= min(abs($candidate1->{exon}->end - $pos_start),abs($candidate1->{exon}->start - $pos_start));
    my $dist2= min(abs($candidate2->{exon}->end - $pos_start),abs($candidate2->{exon}->start - $pos_start));
    if ($dist1 > $dist2) {
      return $candidate2;
    } elsif ($dist1 < $dist2) {
      return $candidate1;
    }
  }
  # If we have not found a better candidate, we use the lexicographic order of the mRNA ID
  my ($mRNA1,$mRNA2) = ($candidate1->{mRNA},$candidate2->{mRNA});
  if(defined $mRNA1 && defined $mRNA1->attribute('ID') && defined $mRNA2 && defined $mRNA2->attribute('ID')) {
    if($mRNA1->attribute('ID') lt $mRNA2->attribute('ID')) {
      return $candidate1;
    } else {
      return $candidate2;
    }
  }
  # If nothing has worked we return "undef"
  return undef;
}




=head1 PRIVATE METHODS

=head2 _init

  Description : init method, load GFF annotation into a
                CracTools::GFF::Query object.

=cut

sub _init {
  my $self = shift;
  my $gff_query;

  # Create a GFF file to query exons
  if($self->mode eq "fast") {
    $gff_query = CracTools::Interval::Query->new();
    my $gff_it = CracTools::Utils::getFileIterator(file => $self->{gff_file},
      parsing_method => sub { CracTools::GFF::Annotation->new(@_) },
      header_regex => "^#",
    );
    while(my $gff_annot = $gff_it->()) {
      $gff_query->addInterval($gff_annot->chr,
        $gff_annot->start+1,
        $gff_annot->end+1,
        $gff_annot->strand,
        $gff_annot,
      );
    }
  } else {
    $gff_query = CracTools::Interval::Query::File->new(file => $self->{gff_file}, type => 'gff');
  }

  $self->{gff_query} = $gff_query;
}

=head2 _constructCandidates

  Arg [1] : String - annot_id
  Arg [2] : Hash ref - candidate
            Since this method is recursive, this is the object that
            we are constructing
  Arg [3] : Hash ref - annot_hash
            annot_hash is a hash reference where keys are annotion IDs
            and values are CracTools::GFF::Annotation objects.

  Description : _constructCandidate is a recursive method that build a
                candidate hash. A candidate is defined as a path into the annotation
                (multi-rooted) tree from a leaf (ex: an exon) to a root (ex: a gene).
  ReturnType  : Candidate Hash ref where keys are GFF features and
                values are CracTools::GFF::Annotation objects :
                { "exon" => CracTools::GFF::Annotation, 
                  "gene" => CracTools::GFF::Annotation,
                  feature => CracTools::GFF::Annotation, ..., 
                  parent_feature => {featureA => featureB},
                  leaf_feature => "exon",
                }

=cut

sub _constructCandidates {
  my ($annot_id,$candidate,$annot_hash) = @_;

  # We init the "leaf_feature" value if this is the first recursion step
  $candidate->{leaf_feature} = $annot_hash->{$annot_id}->feature if !defined $candidate->{leaf_feature};

  my @candidates;
  if (!defined $annot_hash->{$annot_id}){
      carp("Missing feature for $annot_id in the gff file");
  }
  $candidate->{$annot_hash->{$annot_id}->feature} = $annot_hash->{$annot_id};
  my $parents = $annot_hash->{$annot_id}->parents;
  if(@$parents) {
    foreach my $parent (@{$parents}) {
      
      #Test to avoid a deep recursion
      if($parent eq $annot_id) {
  carp("Parent could not be the candidat itself, please check your gff file for $annot_id");
  next;
      # If there is already a parent with this feature type we duplicated
      # the candidate since we are branching in the annotation tree
      }elsif(!defined $annot_hash->{$parent}) {
        carp("Parent not found, please check your gff file for $annot_id (Parent: $parent)");
      
      }elsif(defined $candidate->{$annot_hash->{$parent}->feature}) {
        my %copy_candidate = %{$candidate}; 
        my %copy_parent_feature = %{$candidate->{parent_feature}};
        $copy_candidate{parent_feature} = \%copy_parent_feature; 
        # We register in parent_feature links
        $copy_candidate{parent_feature}->{$annot_hash->{$annot_id}->feature} = $annot_hash->{$parent}->feature;
        my $copy_ref = \%copy_candidate;
        push(@candidates,@{_constructCandidates($parent,$copy_ref,$annot_hash)});
      # If not we only go up to the parent node in order to continue candidate
      # construction
      } else {
        # We register in parent_feature links
        $candidate->{parent_feature}->{$annot_hash->{$annot_id}->feature} = $annot_hash->{$parent}->feature;
        push(@candidates,@{_constructCandidates($parent,$candidate,$annot_hash)});
      }
    }
    return \@candidates;
  } else {
    return [$candidate];
  }
}


=head2 _constructCandidatesFromAnnotation

  Arg [1] : Hash ref - annotations
            Annotions is a hash reference where keys are coordinates
            given by CracTools::Interval::Query::File objects.
  Description : _constructCandidate is a recursive method that build a
                candidate hash.
  ReturnType  : Candidate array ref of all candidates built by _constructCandidate

=cut
sub _constructCandidatesFromAnnotation {
  my $self = shift;
  my $annotations = shift;
  my %annot_hash = ();
  my @candidates = ();

  # Construct annotation hash with annot ID as key
  foreach my $annot_line (@{$annotations}) {
    if($self->mode eq "fast") {
      $annot_hash{$annot_line->attribute('ID')} = $annot_line;
    } else {
      my $annot = CracTools::GFF::Annotation->new($annot_line,'gff3');
      $annot_hash{$annot->attribute('ID')} = $annot;
    }
  }

  # Find leaves in annotation tree
  my %hash_leaves; 
  foreach my $annot_id (keys %annot_hash) {
    #my @parents = $annot_hash{$annot_id}->parents;
    foreach my $parent (@{$annot_hash{$annot_id}->parents}){
      $hash_leaves{$parent} = 1 unless (defined $hash_leaves{$parent});
    }
  }
  foreach my $annot_id (keys %annot_hash) {
    # check if annot_id is a leaf
    if (!defined $hash_leaves{$annot_id}){
      # Get all possible path from this leaf to the root
      push @candidates, @{_constructCandidates($annot_id,my $new_candidate,\%annot_hash)};
    }
  }

  return \@candidates;
}

1;
