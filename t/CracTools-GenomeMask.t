use strict;
use warnings;

use Test::More tests => 15;
use CracTools::GenomeMask;
use CracTools::SAMReader;

use File::Temp 0.23;
use Inline::Files 0.68;  

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

my @pos_set = @{$g_mask->getPosSetInRegion("chr2",2,10)};
is($pos_set[0],4);
is($pos_set[1],8);

is($g_mask->rank("chr2",5),5);
my ($chr,$pos) = $g_mask->select(7);
is($chr,"chr2");
is($pos,8);
is($g_mask->rank("chr2",8),6);

my $sam_file = new File::Temp( SUFFIX => '.sam', UNLINK => 1);
while(<SAM>) {print $sam_file $_;}
close $sam_file;
$g_mask = CracTools::GenomeMask->new(sam_reader => CracTools::SAMReader->new($sam_file));
is($g_mask->getChrLength("2"),242193529,"sam_reader constructor");

my $conf_file = new File::Temp( SUFFIX => '.conf', UNLINK => 1);
while(<CRACCONF>) {print $conf_file $_;}
close $conf_file;
$g_mask = CracTools::GenomeMask->new(crac_index_conf => $conf_file);
is($g_mask->getChrLength("2"),242193529,"crac_index_conf constructor");

__SAM__
@HD	VN:1.5	SO:coordinate
@RG	ID:1	DT:2015-07-13T16:06:12
@PG	ID:1	PN:crac	VN:2.0.0	CL:crac -k 22 --bam -o - --stranded -r reads/ENCSR109IQO_rep1_1.fastq.gz -i /data/indexes/crac/GRCh38 --nb-threads 15 reads/ENCSR109IQO_rep1_2.fastq.gz
@SQ	SN:1	LN:248956422
@SQ	SN:2	LN:242193529
@SQ	SN:3	LN:198295559
@SQ	SN:4	LN:190214555
@SQ	SN:5	LN:181538259
@SQ	SN:6	LN:170805979
@SQ	SN:7	LN:159345973
@SQ	SN:8	LN:145138636
@SQ	SN:9	LN:138394717
@SQ	SN:10	LN:133797422
@SQ	SN:11	LN:135086622
@SQ	SN:12	LN:133275309
@SQ	SN:13	LN:114364328
@SQ	SN:14	LN:107043718
@SQ	SN:15	LN:101991189
@SQ	SN:16	LN:90338345
@SQ	SN:17	LN:83257441
@SQ	SN:18	LN:80373285
@SQ	SN:19	LN:58617616
@SQ	SN:20	LN:64444167
@SQ	SN:21	LN:46709983
@SQ	SN:22	LN:50818468
@SQ	SN:X	LN:156040895
@SQ	SN:Y	LN:57227415
__CRACCONF__
24
1
248956422
2
242193529
3
198295559
4
190214555
5
181538259
6
170805979
7
159345973
8
145138636
9
138394717
10
133797422
11
135086622
12
133275309
13
114364328
14
107043718
15
101991189
16
90338345
17
83257441
18
80373285
19
58617616
20
64444167
21
46709983
22
50818468
X
156040895
Y
57227415
