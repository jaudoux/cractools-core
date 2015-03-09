package CracTools::Interval::Query;
# ABSTRACT: Store and query genomics intervals.
#

use strict;
use warnings;

use CracTools::Utils;
use Set::IntervalTree 0.10;
use Carp;

=head1 METHODS

=head2 new

  Example     : my $intervalQuery = CracTools::Interval::Query->new();
  Description : Create a new CracTools::Interval::Query object
  ReturnType  : CracTools::Interval::Query
  Exceptions  : none

=cut

sub new {
  my $class = shift;

  my %args = @_;

  my $self = bless {
    interval_trees => {},
    }, $class;

  return $self;
}

=head2 addInterval

  Arg [1] : String $chr
            The name of the sequence region that the slice will be
            created on.
  Arg [2] : int $start
            The start of the slice on the sequence region
  Arg [3] : int $end
            The end of the slice on the sequence region
  Arg [4] : int $strand
            The orientation of the slice on the sequence region (1 or -1)
  Arg [5] : Any scalar $value
            The value to be hold by this interval. It can be anything,
            an integer, a hash reference, an array reference, ...

  Example     :
  Description : Add a new genomic interval, with an associated value to the structure

=cut

sub addInterval {
  my $self = shift;
  my ($chr,$start,$end,$strand,$value) = @_;

  my $interval_tree = $self->_getIntervalTree($chr,$strand);
  # If there is no already existing IntervalTree for this ("chr","strand") pair
  if(!defined $interval_tree) {
    # We create a new one
    $interval_tree = Set::IntervalTree->new;
    # We add this new interval tree with the others
    $self->_addIntervalTree($chr,$strand,$interval_tree);
  }

  # We insert the given interval in the IntervalTree
  # pos_end +1 because Interval tree use [a,b) intervals
  #$interval_tree->insert($value,$start,$end+1);
  $interval_tree->insert({value => $value, start => $start, end => $end},$start,$end+1);
}

=head2 fetchByRegion

  Arg [1] : String $seq_region_name
            The name of the sequence region that the slice will be
            created on.
  Arg [2] : int $start
            The start of the slice on the sequence region
  Arg [3] : int $end
            The end of the slice on the sequence region
  Arg [4] : (Optional) int $strand
            The orientation of the slice on the sequence region
  Arg [5] : (Optional) Boolean $windowed
            Only return lines whose intervals are completely contained
            if the slide.

  Example     : my @lines = $IntervalQuery->fetchByRegion('1',298345,309209,'+');
  Description : Retrives lines that belong to the region.
  ReturnType  : Reference to an Array of strings
  Exceptions  : none

=cut

sub fetchByRegion {
  my ($self,$chr,$pos_start,$pos_end,$strand,$windowed) = @_;

  my $interval_tree = $self->_getIntervalTree($chr,$strand);
  
  if(defined $interval_tree) {
    if(defined $windowed && $windowed) {
      # pos_end +1 because Interval tree use [a,b) intervals
      return $self->_processReturnValues($interval_tree->fetch_window($pos_start,$pos_end+1));
    } else {
      # pos_end +1 because Interval tree use [a,b) intervals
      return $self->_processReturnValues($interval_tree->fetch($pos_start,$pos_end+1));
    }
  }
  return [];
}

=head2 fetchByLocation

  Arg [1] : String $seq_region_name
            The name of the sequence region that the slice will be
            created on.
  Arg [2] : int $position
            Location to look for
  Arg [3] : int $strand
            The orientation of the slice on the sequence region

  Example     : my @lines = $intervalQuery->fetchByLocation('1',298345,'+');
  Description : Retrives lines that overlapped the given location.
  ReturnType  : Reference to an Array of strings
  Exceptions  : none

=cut

sub fetchByLocation {
  my ($self,$chr,$position,$strand) = @_;
  return $self->fetchByRegion($chr,$position,$position,$strand);
}

=head2 fetchNearestDown

Search for the closest interval in downstream that does not contain the query
and returns the line associated to this interval. 

=cut

sub fetchNearestDown {
  my ($self,$chr,$position,$strand) = @_;

  my $interval_tree = $self->_getIntervalTree($chr,$strand);
  
  if(defined $interval_tree) {
    my $nearest_down = $interval_tree->fetch_nearest_down($position);
    return ({start => $nearest_down->{start}, end => $nearest_down->{end}},
      $self->_processReturnValue($nearest_down->{value})
    );
  }
  return [];
}

=head2 fetchNearestUp

Search for the closest interval in downstream that does not contain the query
and returns the line associated to this interval. 

=cut

sub fetchNearestUp {
  my ($self,$chr,$position,$strand) = @_;

  my $interval_tree = $self->_getIntervalTree($chr,$strand);
  
  if(defined $interval_tree) {
    my $nearest_up = $interval_tree->fetch_nearest_up($position);
    return ({start => $nearest_up->{start}, end => $nearest_up->{end}},
      $self->_processReturnValue($nearest_up->{value})
    );
  }
  return [];
}

=head2 fetchAllNearestDown

Search for the closest intervals in downstream that does not contain the query
and returns an array reference of lines associated to these intervals. 

=cut

sub fetchAllNearestDown {
  my ($self,$chr,$position,$strand) = @_;

  my ($nearest_down_interval,$nearest_down) = $self->fetchNearestDown($chr,$position,$strand); 
  if(defined $nearest_down) {
    # We return all lines that belong to this
    my ($hits_interval,$hits) = $self->fetchByLocation($chr,$nearest_down_interval->{end},$strand);
    my @valid_hits;
    my @valid_hits_interval;
    for (my $i = 0; $i < @$hits; $i++) {
      # if this inteval as the same "end" boudaries as the nearest down interval
      if($hits_interval->[$i]->{end} == $nearest_down_interval->{end}) {
        push @valid_hits, $hits->[$i];
        push @valid_hits_interval, $hits_interval->[$i];
      }
    }
    return (\@valid_hits_interval,\@valid_hits);
  }
  return [];
}

=head2 fetchAllNearestUp

Search for the closest intervals in downstream that does not contain the query
and returns an array reference of lines associated to these intervals. 

=cut

sub fetchAllNearestUp {
  my ($self,$chr,$position,$strand) = @_;

  my ($nearest_up_interval,$nearest_up) = $self->fetchNearestUp($chr,$position,$strand); 

  if(defined $nearest_up) {
    # We return all lines that belong to this
    my ($hits_interval,$hits) = $self->fetchByLocation($chr,$nearest_up_interval->{start},$strand);
    my @valid_hits;
    my @valid_hits_interval;
    for (my $i = 0; $i < @$hits; $i++) {
      # if this inteval as the same "end" boudaries as the nearest down interval
      if($hits_interval->[$i]->{start} == $nearest_up_interval->{start}) {
        push @valid_hits, $hits->[$i];
        push @valid_hits_interval, $hits_interval->[$i];
      }
    }
    return (\@valid_hits_interval,\@valid_hits);
  }

  return [];
}

=head1 PRIVATE METHODS

=head2 _getIntervalTree 

  $self->_getIntervalTree($chr,$strand);

Return the Set::IntervalTree reference for the chromosome and strand (Default : 1)

=cut

sub _getIntervalTree {
  my ($self,$chr,$strand) = @_;
  $strand = 1 if !defined $strand;
  return $self->{interval_trees}{_getIntervalTreeKey($chr,$strand)};
}

=head2 _addIntervalTree
  
  $self->_addIntervalTree($chr,$strnad,$interval_tree);

Add an Set::IntervalTree object for a specific ("chr","strand") pair.

=cut

sub _addIntervalTree {
  my ($self,$chr,$strand,$interval_tree) = @_;
  $strand = 1 if !defined $strand;
  $self->{interval_trees}{_getIntervalTreeKey($chr,$strand)} = $interval_tree;
}

=head2 _getIntervalTreeKey

  _getIntervalTreeKey($chr,$strand);

Static method that return and unique key for the ("chr","strand") pair passed in arguements.

=cut

sub _getIntervalTreeKey {
  my ($chr,$strand) = @_;
  $strand = 1 if !defined $strand;
  return "$chr"."@"."$strand";
}

=head2 _processReturnValues

  $self->_processReturnValues($array_ref)

Call _processReturnValue() method on each values of the array ref passed in parameters.

=cut

sub _processReturnValues {
  my $self = shift;
  my $return_values = shift;
  my @processed_return_values = ();
  my @processed_return_intervals = ();
  foreach (@{$return_values}) {
    push(@processed_return_values, $self->_processReturnValue($_->{value}));
    push(@processed_return_intervals, { 
        start => $_->{start},
        end => $_->{end}
      }
    );
  }
  return (\@processed_return_intervals,\@processed_return_values);
}

=head2 _processReturnValue

  $self->_processReturnValue($val)

This method process the values contains by each intervals that match a query before returning it.
It is designed to be overloaded by doughter classes.

=cut

sub _processReturnValue {
  my $self = shift;
  my $val = shift;
  return $val;
}

1;
