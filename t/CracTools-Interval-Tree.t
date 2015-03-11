#! /usr/bin/perl

use Test::More tests => 1;
use CracTools::Interval::Tree;
use Set::IntervalTree 0.10;
use Data::Dumper;

use constant NB_LOOP => 1000000;

my $tree = CracTools::Interval::Tree->new();

is(@{$tree->search(1,2)},0);

#for(my $i=0; $i<NB_LOOP; $i++) {
##for(my $i=0; $i<10; $i++) {
#  my $start = int(rand(1000000));
#  my $end = $start + int(rand(1000));
#  $tree->add($start,$end);
#}

#my $tree2 = Set::IntervalTree->new;
#for(my $i=0; $i<NB_LOOP; $i++) {
#  my $start = int(rand(1000000));
#  my $end = $start + int(rand(1000));
#  $tree2->insert($start,$end);
#}
#
#for(my $i=0; $i<100000; $i++) {
#  $tree->search(int(rand(1000000)))
#}
#$tree->add(1,8,'a');
#$tree->add(12,23,'b');
#$tree->add(3,5,'c');
#$tree->add(4,6,'d');
#$tree->add(14,21,'e');
#$tree->add(2,3,'f');
#
#print STDERR Dumper($tree);
#
#my @results = @{$tree->search(2,2)};
#foreach (@results) {
#  print STDERR "$_\n";
#}
