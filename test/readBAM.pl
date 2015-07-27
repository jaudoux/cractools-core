#! /usr/bin/perl
#
use Inline C;
use strict;
use warnings;
#use File::Slurp;

our @cigar_conversion_table = qw(M I D N S H P = X);

my $bam_file = shift;
open(my $fh, "gunzip -c $bam_file |") or die("Cannot open $bam_file");
#
binmode $fh;

my $buffer;

# FIRST READ HEADER

my ($magic,$length) = readBytes($fh,8,$buffer,"A4I");
my ($header_text) = readBytes($fh,$length,$buffer,"a$length");
my ($nb_ref) = readBytes($fh,4,$buffer,"I");

print STDERR "Magic: $magic\n";
print STDERR "Length: $length\n";
print STDERR "# reference: $nb_ref\n";

# Loop over references
for(my $i = 0; $i < $nb_ref; $i++) {
  my ($l_name) = readBytes($fh,4,$buffer,"I");
  my ($name,$l_ref) = readBytes($fh,$l_name+4,$buffer,"Z$l_name I");
  print STDERR "L_name: $l_name, Name: $name, L_ref: $l_ref\n";
}

# NOW WE CAN READ THE RECORDS

my $core_size = 32;
while(my ($block_size) = readBytes($fh,4,$buffer,"I")) {
  my $remain_size = $block_size - $core_size;
  read($fh,$buffer,$block_size);

  my($ref_id,$pos,$bin_mq_nl,$flag_nc,$l_seq,$next_ref_id,$next_pos,$tlen,$remain)=unpack("IIVVIIIIa$remain_size",$buffer);

  my $bin = $bin_mq_nl >> 16;
  my $mapq = $bin_mq_nl >> 8 & 0xff;
  my $l_read_name = $bin_mq_nl & 0xff;
  my $flag = $flag_nc >> 16;
  my $n_cigar = $flag_nc & 0xffff;

  #my ($read_name) = readBytes($fh,$l_read_name,$buffer,"Z$l_read_name");
  #
  my $seq_byte_length = ($l_seq+1)*0.5;
  $remain_size -= ($l_read_name + 4*$n_cigar + $seq_byte_length + $l_seq);
  #my ($read_name,$cigar,$seq,$qual,$remain) = unpack("Z$l_read_name V$n_cigar C$seq_byte_length A$l_seq a$remain_size",$remain);
  my ($read_name,@values) = unpack("Z$l_read_name V$n_cigar C$seq_byte_length A$l_seq a$remain_size",$remain);
  #my @cigar = splice (@values,0,$n_cigar);
  #unpackCigar(\@cigar,$n_cigar);
  #my @seq_bytes = splice (@values,0,$seq_byte_length);
  #my $qual = shift @values;
  print STDERR "$read_name\t\n";#.unpackCigar(\@cigar,$n_cigar)."\n";
  sleep 1;
}

sub readBytes {
  my ($fh,$nb_char,$buffer,$template) = @_;
  read($fh,$buffer,$nb_char);
  return unpack($template,$buffer);
}

sub unpackCigar {
  my ($packed_cigar, $cigar_length) = @_;
  my $cigar;
  for(my $i = 0; $i < $cigar_length; $i++) {
    my $l = $packed_cigar->[$i] >> 4;
    my $op = $cigar_conversion_table[$packed_cigar->[$i] & 0xf];
    $cigar .= $l.$op;
  }
  return $cigar;
}

__END__
__C__

