package CracTools::GenomeMask;
# ABSTRACT: A bit vector mask over the whole genome
#

use strict;
use warnings;

use CracTools::BitVector;
use CracTools::Utils;
use Carp;

=head1 METHODS

=head2 new

There is mutiple ways to create a genome mask:

One can specify a argument called C<genome> that is a hashref where keys are chromosome names
and values are chromosomes length.

  my $genome_mask = CracTools::GenomeMask->new( genome => { seq_name => length,
                                                            seq_name => length,
                                                            ...} );
One can specify a argument called C<crac_index_conf> that the configuration file of a CRAC index

  my $genome_mask = CracTools::GenomeMask->new(crac_index_conf => file.conf);

One can specify a C<CracTools::SAMReader> object in order to read chromosomes names and lenght from
the header

  my $genome_mask = CracTools::GenomeMask->new(sam_reader => CracTools::SAMReader->new(file.sam));


=cut

sub new {
  my $class = shift;

  my %args = @_;

  # the genome mask is not stranded by default
  my $is_stranded = defined $args{is_stranded}? $args{is_stranded} : 0;
  my $verbose = defined $args{verbose}? $args{verbose} : 0;

  my %bit_vectors;

  if(defined $args{genome}) {

    foreach my $chr (keys %{$args{genome}}) {
      if(!defined $bit_vectors{$chr}) {
        $bit_vectors{$chr} = CracTools::BitVector->new($args{genome}->{$chr});
      } else {
        croak "Multiple definition of sequence $chr inf the genome lengths";
      }
    }
  } elsif(defined $args{crac_index_conf}) {
    print STDERR "Creating GenomeMask from crac index conf file : $args{crac_index_conf}\n" if $verbose;
    my $conf_fh = CracTools::Utils::getReadingFileHandle($args{crac_index_conf});
    my $nb_chr = <$conf_fh>;
    my $nb_chr_found = 0;
    while(<$conf_fh>) {
      my $chr = $_;
      chomp $chr;
      my $chr_length = <$conf_fh>;
      chomp $chr_length;
      if(defined $chr_length) {
        print STDERR "\tCreating bitvecor for chr $chr of length $chr_length\n" if $verbose;
        $bit_vectors{$chr} = CracTools::BitVector->new($chr_length);
        $nb_chr_found++;
      } else {
        croak "Missing genome length for chromosome $chr";
      }
    }
    croak "There is less chromosome found ($nb_chr_found) in $args{crac_index_conf} than expected ($nb_chr)" if $nb_chr_found < $nb_chr;
  } elsif(defined $args{sam_reader}) {
    my $refseq_lengths = $args{sam_reader}->allRefSeqLengths();
    foreach my $chr (keys %{$refseq_lengths}) {
      $bit_vectors{$chr} = CracTools::BitVector->new($refseq_lengths->{$chr});
    }
  } else {
    croak "There is no valid argument to extract the chromosomes names and length";
  }

  my $self = bless {
    bit_vectors => \%bit_vectors,
    is_stranded => $is_stranded,
  }, $class;

  return $self;
}

sub isStranded {
  my $self = shift;
  return $self->{is_stranded};
}

=head2 getBitvector

Return the C<CracTools::BitVector> associated with the reference name given in argument.

If no bitvectors exists for this reference, a warning will be reported.

=cut

sub getBitvector {
  my $self = shift;
  my $chr = shift;
  my $strand = shift;
  if(defined $self->{bit_vectors}->{$chr}) {
    return $self->{bit_vectors}->{$chr};
  } else {
    carp "There is no bitvector for sequence $chr in the genome mask";
    return undef;
  }
}

=head2 getChrLength

=cut

sub getChrLength {
  my $self = shift;
  my $chr = shift;
  my $bv = $self->getBitvector($chr);
  return defined $bv? $bv->length : undef;
}

=head2 setPos

  $genome_mask->setPos($chr,$pos)

Set the bit to 1 for this position

=cut

sub setPos {
  my ($self,$chr,$pos) = @_;
  my $bv = $self->getBitvector($chr);
  $bv->set($pos) if defined $bv;
}

=head2 setRegion

  $genome_mask->setRegion($chr,$start,$end)

Set all bits to 1 for this region

=cut

sub setRegion {
  my ($self,$chr,$start,$end) = @_;
  for(my $i = $start; $i <= $end; $i++) {
    $self->setPos($chr,$i);  
  }
}

=head2 getPos

  my $boolean = $genome_mask->getPos($chr,$pos)

Retrun true is the bit is set for this position

=cut
 
sub getPos {
  my ($self,$chr,$pos) = @_;
  my $bv = $self->getBitvector($chr);
  return $bv->get($pos) if defined $bv;
}

=head2 getPosSetInRegion

  my $nb_pos_set = $genome_mask->getNbBitsSetInRegion($chr,$start,$end)

Return the positions of bits set in this genomic region

=cut

sub getPosSetInRegion {
  my ($self,$chr,$start,$end) = @_;
  my $bv = $self->getBitvector($chr);
  my @pos;
  if(defined $bv) {
    for(my $i = $start; $i <= $end; $i++) {
      push(@pos,$i) if $bv->get($i) == 1;  
    }
  }
  return \@pos;
}

=head2 getNbBitsSetInRegion

  my $nb_pos_set = $genome_mask->getNbBitsSetInRegion($chr,$start,$end)

Return the number of bit set in this genomic region

=cut

sub getNbBitsSetInRegion {
  my $self = shift;
  return scalar @{$self->getPosSetInRegion(@_)};
}

=head2 rank

Return the number of bit set in the genome before this position

=cut

sub rank {
  my ($self,$chr,$pos) = @_;
  my $cumulated_bits = 0;
  my $i = 0;
  my @chr_sorted = sort keys %{$self->{bit_vectors}};
  while($chr_sorted[$i] ne $chr) {
    $cumulated_bits += $self->getBitvector($chr_sorted[$i])->nb_set;
    $i++;
  }
  return $cumulated_bits + $self->getBitvector($chr_sorted[$i])->rank($pos);
}

=head2 select 

Return an array with a (chr,pos) of the Nth bit set

=cut

sub select {
  my $self = shift;
  my $i = shift;
  my $cumulated_bits = 0;
  my @chr_sorted = sort keys %{$self->{bit_vectors}};
  my $j = 0;
  while($j < @chr_sorted && $cumulated_bits + $self->getBitvector($chr_sorted[$j])->nb_set < $i) {
    my $chr = $chr_sorted[$j];
    my $bv = $self->getBitvector($chr);
    $cumulated_bits += $self->getBitvector($chr)->nb_set;
    $j++;
  }
  my $chr = $chr_sorted[$j-1];
  my $pos = $self->getBitvector($chr_sorted[$j-1])->select($i - $cumulated_bits + 1);
  return ($chr,$pos);
}




1;
