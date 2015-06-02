#! /usr/bin/perl
#

use strict;
use warnings;
use CracTools::Aligner::Crac;
use Data::Dumper;

#Fasta file
my @fasta= "../extra/toy_1.fasta ../extra/toy_2.fasta";

my $obj = CracTools::Aligner::Crac->new(@fasta);

print Dumper($obj);
