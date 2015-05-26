#! /usr/bin/perl

use Test::More tests => 11;
use CracTools::BitVector;

my $bv = CracTools::BitVector->new(10);

# Length
is($bv->length(),10);

# Nb_set
is($bv->nb_set(),0);

$bv->set(2);
$bv->set(5);

#0010010100


is($bv->nb_set(),2);

# get/set
is($bv->get(2),1);
is($bv->get(1),0);

# unset
$bv->unset(5);
is($bv->get(5),0);

# prev
is($bv->prev(5),2);

# succ
$bv->set(7);
is($bv->succ(3),7);

#0010000100

# rank
is($bv->rank(5),1);
is($bv->rank(7),2);

# select
is($bv->select(2),7);
