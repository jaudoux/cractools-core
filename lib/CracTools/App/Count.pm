package CracTools::App::Count;
# ABSTRACT: Count reads based on a GFF/GTF annotation file

use strict;
use warnings;

use Carp;
use CracTools;
use CracTools::Utils; # For multi-purpose tasks
use CracTools::Annotator; # This help to index the annotation file for fast query
use Data::Dumper; # Debugging purpose

sub new {
  my $class = shift;

  my %args = @_;

  my $annot_file = $args{gff_file};

  croak "Missing gff_file argument" unless defined $annot_file;

  # Create an annotator object to query the GFF annotation file
  my $annotator = CracTools::Annotator->new($annot_file,"fast");

  my $self = bless {
    feature_type => $args{feature_type},
    is_stranded => $args{is_stranded},
    annotator => $annotator,
  }, $class;

  return $self;
}

sub isStranded {
  my $self = shift;
  return $self->{is_stranded};
}

sub featureType {
  my $self = shift;
  return $self->{feature_type};
}

sub annotator {
  my $self = shift;
  return $self->{annotator};
}

sub getCounts {
  my $self = shift;
  my $bam_file = shift;

  my $annotator = $self->annotator;

  my $bam_it = CracTools::Utils::bamFileIterator($bam_file,"");

  my %counts;
  my $prev_pos = 0;
  my $prev_cigar = "";
  my @prev_candidates = ();

  while(my $line = $bam_it->()) {
    # We parse the sam line with a lightweigh method
    my $sam_line = CracTools::Utils::parseSAMLineLite($line);
    if($prev_pos == $sam_line->{pos} && $prev_cigar eq $sam_line->{original_cigar}) {
      foreach my $cand (@prev_candidates) {
        $counts{$cand}++;
      }
      # We go to next sam line and avoid the counting process
      next;
    } else {
      @prev_candidates = ();
      $prev_pos = $sam_line->{pos};
      $prev_cigar = $sam_line->{original_cigar};
    }
    #print STDERR "Read: ".$sam_line->{qname}."\n";

    # #################
    # 1. Find read chunks and register for each on of them a list of candidates
    # -1 because annotator use a 0-basede coordinate system
    my $low = $sam_line->{pos} - 1;
    my $high = $low;
    my $strand = $sam_line->{flag} & 16? -1 : 1;
    # This table will hold the candidates for each chunk of the read
    my @candidates;
    my $nb_chunk = 0;
    foreach my $cigar_chunk (@{$sam_line->{cigar}}) {
      # If we have a N operator in the sequence it means that we have a splited
      # read that may correspond to a splice junction
      if($cigar_chunk->{op} =~ 'N') {
        # If this is the first chunk we see, we try to find candidates
        # that have exon's end close to our chunk's end
        if($nb_chunk == 0) {
          # This chunk's end should be close to an exon end
          push(@candidates,$self->getReadChunkCandidates($annotator,
              $sam_line->{rname},
              $low,
              $high,
              $strand,
              \&compareSubExonEnd
            )
          );
          #print STDERR Dumper(\@candidates);
        } else {
          # If this is not the first chunk, then we have a middle chunk
          # that should correspond to a whole exon
          # This chunk's boundaries should correspond to an exon
          push(@candidates,$self->getReadChunkCandidates($annotator,
              $sam_line->{rname},
              $low,
              $high,
              $strand,
              \&compareSubExon
            )
          );
        }
        # We adjust boundaries to start the new chunk
        $low = $high + $cigar_chunk->{nb};
        $high = $low;
        $nb_chunk++;
      } elsif($cigar_chunk->{op} !~ /[SHI]/) {
        # We move the upper bound further
        #print STDERR "cigar op: ".$cigar_chunk->{op}."\n";
        $high = $high + $cigar_chunk->{nb};
      }
    }
    # We have reach the cigar, we can count the last chunk
    # If this read is not spliced, we try to find a candidate where
    # the read ins included in an exon
    if($nb_chunk == 0) {
      push(@candidates,$self->getReadChunkCandidates($annotator,
          $sam_line->{rname},
          $low,
          $high,
          $strand,
          \&compareSubExonIncluded
        )
      );
    } else {
      # This chunk is a last one of a spliced read, it should
      # correspond to an exon start
      push(@candidates,$self->getReadChunkCandidates($annotator,
          $sam_line->{rname},
          $low,
          $high,
          $strand,
          \&compareSubExonStart)
      );
    }
    
    # #################
    # 2. We loop over the candidates by comparing candidates of adjacent
    #    chunks and finding the best possible candidates:
    #
    #    - If two adjacent chunk contains candidate with the same feature
    #      that we are conting
    #    - If there is no such feature we look for the closest mutual parent
    #      of each candidate, and keeps only
    my $min_dist;
    # First we init the hash of selected candidates for the first chunk
    my %counted_candidates;
    my %selected_candidates;
    if(@candidates < 2) {
      foreach my $first_chunk_cand (@{$candidates[0]}) {
        $selected_candidates{$first_chunk_cand->{$self->featureType}->attribute('ID')} = $first_chunk_cand->{$self->featureType};
      }
    }
    # We loop over all chunk candidates to find the best annotation feater
    # by selecting those that have the same exons.
    for(my $i = 1; $i < @candidates; $i++) {
      my %current_selection;
      my $current_min_dist;# = $min_dist;
      # We do the intersection with the previous chunk
      foreach my $prev_cand (@{$candidates[$i-1]}) {
        foreach my $curr_cand (@{$candidates[$i]}) {
          my $prev_parent_feat = $self->featureType;
          my $curr_parent_feat = $self->featureType;

          my $dist = 0;
          # We declare a hash that stores "already seen features" in order
          # to avoid infinite loop.
          my %avoid_recurs;
          # We loop over parent feature to find the distance from the closest parent
          while(defined $prev_parent_feat && defined $curr_parent_feat &&
            $prev_cand->{$prev_parent_feat}->attribute('ID') ne $curr_cand->{$curr_parent_feat}->attribute('ID')) {
            $dist++;
            if(defined $avoid_recurs{$prev_parent_feat} || 
              defined $avoid_recurs{$curr_parent_feat}) {
              carp "GFF deep recursion\n";
              last;
            }
            $avoid_recurs{$prev_parent_feat} = 1;
            $avoid_recurs{$curr_parent_feat} = 1;
            $prev_parent_feat = $prev_cand->{parent_feature}->{$prev_parent_feat};
            $curr_parent_feat = $curr_cand->{parent_feature}->{$curr_parent_feat};
          }
          if(!defined $current_min_dist || $dist < $current_min_dist) {
            $current_min_dist = $dist;
            %current_selection = ();
          }
          if($dist == $current_min_dist) {
            $current_selection{$curr_cand->{$self->featureType}->attribute('ID')} = $curr_cand->{$self->featureType};
            $current_selection{$prev_cand->{$self->featureType}->attribute('ID')} = $prev_cand->{$self->featureType};
          }
        }
      }
      # If we have something in the current selection; that means that we have
      # found a good candidate(s)
      if(%current_selection && (!defined $min_dist || $current_min_dist == $min_dist)) {
        # We merge candidate hashes
        @selected_candidates{keys %current_selection} = values %current_selection;
        $min_dist = $current_min_dist;
      # If we do not, we need to count the read for the prev candidate(s)
      # and continue for next chunks to come
      } else {
        foreach my $candidate (keys %selected_candidates) {
          $counted_candidates{$candidate} = 1;
        }
        %selected_candidates = %current_selection;
        $min_dist = $current_min_dist;
      }
    }
    # We add to the count table what's left
    foreach my $candidate (keys %selected_candidates) {
      $counted_candidates{$candidate} = 1;
    }

    # Now we count all candidates the we have encounter
    foreach my $candidate (keys %counted_candidates) {
      $counts{$candidate}++;
      push(@prev_candidates,$candidate);
    }
  }

  return \%counts;
}

sub getReadChunkCandidates {
  my $self = shift;
  my ($annotator,$chr,$start,$end,$strand,$compareSub) = @_;
  my @potential_candidates;
  my @candidates;
  my $is_included = 0;
  if($self->isStranded && defined $strand) {
    # TODO
  } else {
    foreach my $query_strand ((1,-1)) {
      my ($best_candidates) = $annotator->getBestAnnotationCandidates($chr,
          $start,
          $end,
          $query_strand,
          \&prioritySub, 
          $compareSub,
      );
      push @potential_candidates, @{$best_candidates};
    }
  }
  foreach my $candidate (@potential_candidates) {
    # First we want the candidate to have the feature type that we are
    # counting on
    if(defined $candidate->{$self->featureType}) {
      # if candidate is included
      if($candidate->{$self->featureType}->start >= $start && $candidate->{$self->featureType}->end <= $end) {
        if(!$is_included) {
          @candidates = ();
          $is_included = 1;
        }
        push(@candidates,$candidate);
      } elsif(!$is_included) {
        push(@candidates,$candidate);
      }
    }
  }
  return \@candidates;
}

sub prioritySub {
  my ($pos_start,$pos_end,$candidate) = @_;
  my ($priority,$type) = (-1,'');
  my $exon = $candidate->{exon};
  if (defined $exon) {
    $priority = 1;
    $type = 'EXON';
  } else {
    $priority = 2;
    $type = 'NO_EXON';
  }
  return ($priority,$type);
}

sub compareSubExonIncluded {
  my ($candidate1,$candidate2,$pos_start,$pos_end) = @_;
  # If both candidates are exons we try to find wich one is closer to the
  # pos_start of the region to annotate
  if ($candidate1->{exon} && $candidate2->{exon}) { 
    my $candidate1_included = $pos_start >= $candidate1->{exon}->start && 
                              $pos_end   <= $candidate1->{exon}->end;
    my $candidate2_included = $pos_start >= $candidate2->{exon}->start &&
                              $pos_end   <= $candidate2->{exon}->end;
    if($candidate1_included && !$candidate2_included) {
      return $candidate1;
    } elsif($candidate2_included && !$candidate1_included) {
      return $candidate2;
    }
  }
  # If nothing has worked we return "undef"
  return undef;
}

# Return the candidate that is the closest to exon bounds
sub compareSubExon {
  my ($candidate1,$candidate2,$pos_start,$pos_end) = @_;
  # If both candidates are exons we try to find wich one is closer to the
  # pos_start of the region to annotate
  if ($candidate1->{exon} && $candidate2->{exon}) { 
    my $dist1 = abs($candidate1->{exon}->start - $pos_start) + abs($candidate1->{exon}->end - $pos_end);
    my $dist2 = abs($candidate2->{exon}->start - $pos_start) + abs($candidate2->{exon}->end - $pos_end);
    if($dist1 < $dist2) {
      return $candidate1;
    } elsif($dist1 > $dist2) {
      return $candidate2;
    }
  }
  # If nothing has worked we return "undef"
  return undef;
}

# Return the candidate that is the closest to exon start bound
sub compareSubExonStart {
  my ($candidate1,$candidate2,$pos_start,$pos_end) = @_;
  # If both candidates are exons we try to find wich one is closer to the
  # pos_start of the region to annotate
  if ($candidate1->{exon} && $candidate2->{exon}) { 
    my $dist1 = abs($candidate1->{exon}->start - $pos_start);
    my $dist2 = abs($candidate2->{exon}->start - $pos_start);
    if($dist1 < $dist2) {
      return $candidate1;
    } elsif($dist1 > $dist2) {
      return $candidate2;
    }
  }
  # If nothing has worked we return "undef"
  return undef;
}

# Return the candidate that is the closest to exon end bound
sub compareSubExonEnd {
  my ($candidate1,$candidate2,$pos_start,$pos_end) = @_;
  # If both candidates are exons we try to find wich one is closer to the
  # pos_end of the region to annotate
  if ($candidate1->{exon} && $candidate2->{exon}) { 
    my $dist1 = abs($candidate1->{exon}->end - $pos_end);
    my $dist2 = abs($candidate2->{exon}->end - $pos_end);
    if($dist1 < $dist2) {
      return $candidate1;
    } elsif($dist1 > $dist2) {
      return $candidate2;
    }
  }
  # If nothing has worked we return "undef"
  return undef;
}

1;
