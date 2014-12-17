use strict;
use warnings;

use Test::More tests => 7;
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
