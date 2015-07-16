package CracTools::BitVector;
# ABSTRACT: Full Perl BitVector implementation

use Exporter;
our @ISA = qw(Exporter);

use strict;
use POSIX qw(floor ceil);
use constant BITNESS => 64;

# use overload 
#     q{=} => 'copy';
#     q{""} => 'to_string',
#     '*' => 'mult',
#     '/' => 'divide',
#     '+' => 'add',
#     '.' => 'multElemPerElem',
#     '-' => 'minus';

# BitVector Constructor
# 1 parameter: length of the vector
sub new {
    my ($class, $length)  = @_;
    my $this = {};
    bless($this, $class);
    
    $this->{n} = $length;
    $this->{set} = 0;
    $this->{first_bit_set} = $length + 1;

    my $max = int($length / BITNESS);
    if ($length % BITNESS > 0) {
        $max++;
    }
    $this->{maxCases} = $max;

    for (my $i=0; $i < $max; $i++) {
        push @{$this->{bits}}, 0;
    }

    return $this;
}

sub firstBitSet {
  my $self = shift;
  return $self->{first_bit_set};
}

sub copy {
    my ($this) = @_;
    my $new = {};
    $new->{n} = $this->{n};
    $new->{set} = $this->{set};
    $new->{maxCases} = $this->{maxCases};
    @{$new->{bits}} = @{$this->{bits}};
    if (defined($this->{block})) {
	@{$new->{block}} = @{$this->{block}};
    }
    bless $new, ref($this);
    return $new;
}

# Set a 1-bit at position i 
# @param i: position in the bit vector
# @post get($i) == 1
sub set {
    my ($this, $i) = @_;
    $this->{first_bit_set} = $i if $i < $this->{first_bit_set};
    if (! $this->get($i)) {
        $this->{set}++;
        $this->{bits}[int($i / BITNESS)] |= (1 << ($i % BITNESS));
        $this->{block}[int($i / BITNESS)]++;
    }
}

# Unset a 1-bit at position i 
# @param i: position in the bit vector
# @post get($i) == 0
sub unset {
    my ($this, $i) = @_;
    if ($this->get($i)) {
        $this->{set}--;
        $this->{bits}[int($i / BITNESS)] &= (~(1 << ($i % BITNESS)));
        $this->{block}[int($i / BITNESS)]--;
    }
}

# @param i position in the bit vector
# @return the value of the bit at position i
sub get {
    my ($this, $i) = @_;
    if (! defined($this) || ! defined($this->{bits})
	|| ! defined($this->{bits}[int($i / BITNESS)])) {
	die("Bad position i=$i case -> ".int($i/BITNESS)."\n");
    }
    return ($this->{bits}[int($i / BITNESS)] & (1 << ($i % BITNESS)))? 1 : 0;
}

# @param i position in the bit vector
# @param max (optional) maximal shift from i
# @return The previous position that has a bit set.
#         ie. a position j <= i such that get(j) == 1
#         && get(k) == 0 , j < k <= i
#         If max is set, j >= i - max
#         -1 if no such position exists
sub prev {
    my ($this, $i, $max) = @_;
    my $start = $i;
    
    # Is there a 1 in the current block?
    if ($this->{bits}[int($i / BITNESS)] 
        & ((~0) >> (BITNESS-1) - ($i % BITNESS))) {
        while ($this->get($i) == 0) {
            $i--;
        }
	if (defined($max) && $i < $start - $max) {
	    return -1;
	}
        return $i;
    } else {
        # Look for a block with a 1-bit set.

        my $block = int($i / BITNESS) - 1;
        while ($block >= 0
               && (! defined($this->{block}[$block])
                   || ! $this->{block}[$block])
	       && (! defined($max) 
		   || $i - ($block+1) * (BITNESS) + 1  <= $max)) {
            $block--;
        }
        if ($block < 0 || (defined($max) 
			   && $i - ($block+1) * (BITNESS) + 1 > $max)) {
            return -1;
        }
        $i = ($block + 1) * (BITNESS) - 1;
        while ($this->get($i) == 0) {
            $i--;
        }
	if (defined($max) && $i < $start - $max) {
	    return -1;
	}	
        return $i;
    }
}

# @param i position in the bit vector
# @param max (optional) maximal shift from i
# @return The next position that has a bit set.
#         ie. a position j >= i such that get(j) == 1
#         && get(k) == 0 , j > k >= i
#         If max is set, j <= i + max
#         -1 if no such position exists
sub succ {
    my ($this, $i, $max) = @_;
    my $start = $i;

    # Is there a 1 in the current block?
    if ($this->{bits}[int($i / BITNESS)] 
        & ((~0) << ($i % BITNESS))) {
        while ($i < $this->{n} && $this->get($i) == 0) {
            $i++;
        }
	if ($i == $this->{n} || (defined($max) && $i > $max + $start)) {
	    return -1;
	}
        return $i;
    } else {
        # Look for a block with a 1-bit set.

        my $block = int($i / BITNESS) + 1;
        while ($block < $this->{maxCases} 
               && (! defined($this->{block}[$block])
                   || ! $this->{block}[$block])
	       && (! defined($max) 
		   || $block * (BITNESS) - $i <= $max)) {
            $block++;
        }
        if ($block >= $this->{maxCases} 
	    || (defined($max) 
		&& $block * (BITNESS) - $i > $max)) {
            return -1;
        }
        $i = $block * (BITNESS);
        while ($i < $this->{n} && $this->get($i) == 0) {
            $i++;
        }
	if ($i == $this->{n} || defined($max) && $i > $max + $start) {
	    return -1;
	}	
        return $i;
    }
}

# @param i position in the bit vector
# @return The next position that has a bit set.
#         ie. a position j >= i such that get(j) == 1
#         && get(k) == 0 , j > k >= i
#         If max is set, j <= i + max
#         -1 if no such position exists
sub rank {
  my ($self, $i) = @_;
  my $rank = $self->get($i)? 1 : 0;
  my $found_bit = 1;
  while($found_bit) {
    $i = $self->prev($i-1);
    if ($i != -1) {
      $rank++;
    } else {
      $found_bit = 0;
    }
  }
  return $rank;
}

# @param i position in the bit vector
# @return The next position that has a bit set.
#         ie. a position j >= i such that get(j) == 1
#         && get(k) == 0 , j > k >= i
#         If max is set, j <= i + max
#         -1 if no such position exists
sub select {
  my ($self, $nb) = @_;
  my $i = $self->firstBitSet;
  while ($nb > 1 && $i != -1) {
    $i = $self->succ($i+1);
    $nb--;
  }
  return $i;
}

sub length {
    my ($this) = @_;
    return $this->{n};
}

sub nb_set {
    my ($this) = @_;
    return $this->{set};
}

sub to_string {
    my $this = shift;
    my $output = '';
    
    for (my $i=0; $i < $this->{n}; $i++) {

        if ($i % BITNESS == 0 && $i > 0) {
            $output .= ' ';
        }
        $output .= $this->get($i);
    }
    return $output
}

1;
