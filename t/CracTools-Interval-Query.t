#! /usr/bin/perl
#
use strict;
use warnings;

use Test::More tests => 17;
use CracTools::Interval::Query;
use File::Temp;
use Inline::Files;

my $gff_line = "1\tEnsembl\texon\t1\t10\t.\t+\t0\n";
my $sam_line = "HWI-ST170:310:8:1101:1477-2249/2\t161\t12\t49333532\t254\t30M217N35M\t12\t49334794\t0\tCGATCATTGCTGTCGACCACAAATATCAACCCTTGGGTGTTCTGGAAGTAGTGTCTCCAGAGGGG\t__[ceccececggfhhfhf^ggagghhiifhifhhf^^eeeeefdfhcf`aabfcddgggf^";
my $short_bed_line = "chr22\t1000\t5000";
my $bed_line = "chr22\t1000\t5000\tcloneA\t960\t+\t1000\t5000\t0\t2\t567,488,\t0,3512";

my $interval = CracTools::Interval::Query::_getIntervalsFromGFFLine($gff_line);

is($interval->[0]{low},1);
is($interval->[0]{high},10);

$interval = CracTools::Interval::Query::_getIntervalsFromSAMLine($sam_line);
is($interval->[0]{low},49333532);
is($interval->[0]{high},49333562);
is($interval->[1]{low},49333779);
is($interval->[1]{high},49333814);

$interval = CracTools::Interval::Query::_getIntervalsFromBEDLine($short_bed_line);
is($interval->[0]{low},1001);
is($interval->[0]{high},5000);

$interval = CracTools::Interval::Query::_getIntervalsFromBEDLine($bed_line);
is($interval->[0]{low},1001);
is($interval->[0]{high},1567);
is($interval->[1]{low},4513);
is($interval->[1]{high},5000);

# Create a temp file with the SAM lines described above
my $gff_file = new File::Temp( SUFFIX => '.gff', UNLINK => 1);
while(<GFF>) {print $gff_file $_;}
#print $gff_file $gff;
close $gff_file;

my $intervalQuery = CracTools::Interval::Query->new(file => $gff_file,
                                               type => 'gff',
                                             );

is(@{$intervalQuery->fetchByLocation(1,3,1)}, 3, 'findAnnotations() (1)');
is(@{$intervalQuery->fetchByLocation(1,3,'-1')}, 0, 'findAnnotations() (2)');
is(@{$intervalQuery->fetchByLocation(1,4,'1')}, 2, 'findAnnotations() (3)');
ok($intervalQuery->fetchByLocation(1,10,'1')->[0] =~ /1\t10\t\.\t\+/, 'findAnnotations() (4)');
is($intervalQuery->fetchByRegion(1,3,3,1)->[0],$intervalQuery->fetchByLocation(1,3,1)->[0],'fetchByRegion');

__GFF__
1	Ensembl	exon	1	10	.	+	0
1	Ensembl	exon	2	3	.	+	0
1	Ensembl	exon	3	9	.	+	0
1	Ensembl	exon	1	2	.	-	0
2	Ensembl	exon	4	8	.	+	0
