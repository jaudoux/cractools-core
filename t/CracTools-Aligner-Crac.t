#! /usr/bin/perl
#

use strict;
use warnings;
use CracTools::Aligner::Crac;
use Data::Dumper;

#Fasta file
my $f = "-r ../extra/toy_1.fasta ../extra/toy_2.fasta";
my $index = "-i /data/indexes/crac/GRCh38";
my $k = "-k 20";
my $options = "--detailed-sam";

my @list = ($f,$index,$k,$options);

my $obj = CracTools::Aligner::Crac->new(@list);

print Dumper($obj);
