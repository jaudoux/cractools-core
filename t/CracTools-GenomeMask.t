use strict;
use warnings;

use Test::More tests => 11;
use CracTools::GenomeMask;

my $g_mask = CracTools::GenomeMask->new(genome => {chr1 => 10, chr2 => 20});

is($g_mask->getChrLength("chr1"),10);
$g_mask->setRegion("chr1",2,5);

is($g_mask->getPos("chr1",1),0);
is($g_mask->getPos("chr1",2),1);
is($g_mask->getPos("chr1",3),1);
is($g_mask->getPos("chr1",5),1);
is($g_mask->getPos("chr1",6),0);

is($g_mask->getNbBitsSetInRegion("chr1",3,6), 3);

$g_mask->setPos("chr2",4);
$g_mask->setPos("chr2",8);

is($g_mask->rank("chr2",5),5);
my ($chr,$pos) = $g_mask->select(7);
is($chr,"chr2");
is($pos,8);
is($g_mask->rank("chr2",8),6);
