use strict;
use warnings;
use Test::More tests => 11;
use CracTools::App::Count;
use CracTools::Annotator;
use Inline::Files 0.68;
use File::Temp;

# Testing GFF1
{
  # Get GFF as temp file
  my $gff_file = new File::Temp( SUFFIX => '.gff', UNLINK => 1);
  while(<GFF>) {print $gff_file $_;}
  close $gff_file;

  # get SAM as temp file
  my $sam_file = new File::Temp( SUFFIX => '.sam', UNLINK => 1);
  while(<SAM>) {print $sam_file $_;}
  close $sam_file;

  my $counter_exon = CracTools::App::Count->new(gff_file => $gff_file, feature_type => "exon", is_stranded => 0);

  my %exons_counts = %{$counter_exon->getCounts($sam_file)};
  #my $counts = CracTools::App::Count::getCounts($gff_file,$sam_file,"exon",0);
  is($exons_counts{E1},3);
  is($exons_counts{E1b},3);
  is($exons_counts{E3},2);
  is($exons_counts{E3b},1);
  is($exons_counts{E2},1);
  is($exons_counts{E4},1);
  is($exons_counts{E5},1);

  my $counter_transcripts = CracTools::App::Count->new(gff_file => $gff_file, feature_type => "mRNA", is_stranded => 0);
  my %transcripts_counts = %{$counter_transcripts->getCounts($sam_file)};
  #my $counts = CracTools::App::Count::getCounts($gff_file,$sam_file,"exon",0);
  is($transcripts_counts{T1},2);
  is($transcripts_counts{T2},3);
  is($transcripts_counts{T3},2);
  is($transcripts_counts{T4},1);
}

# GFF1 structure
#
# Gene G1:
#               111111111122222222223333333
#      123456789012345678901234567890123456
#               +----------------------+
#               |          G2          |
#               +----------------------+
# T4 :          +---------+      +----+
#               |    E4   |------| E5 |
#               +---------+      +----+
#      +----------------------------------+
#      |               G1                 |
#      +----------------------------------+
# T1 : +----+                        +----+
#      | E1 |------------------------| E3 |
#      +----+                        +----+
# T2 : +-------+                     +----+
#      | E1b   |---------------------| E3 |
#      +-------+                     +----+
# T3 : +----+        +----+        +------+
#      | E1 |--------| E2 |--------|  E3b |
#      +----+        +----+        +------+
# read1  ###                              
# read2   ####                            
# read3     ####---------------------###
# read4  ####--------#####---------###
# read5  ####-------------------------###
# read6               #####------####


__GFF__
1	line1	exon	1	6	.	+	0	ID=E1;Parent=T1,T3
1	line2	exon	1	9	.	+	0	ID=E1b;Parent=T2
1	line1	exon	14	20	.	+	0	ID=E2;Parent=T3
1	line1	exon	30	36	.	+	0	ID=E3;Parent=T1,T2
1	line1	exon	28	36	.	+	0	ID=E3b;Parent=T3
1	line1	mRNA	1	36	.	+	0	ID=T1;Parent=G1
1	line1	mRNA	1	36	.	+	0	ID=T2;Parent=G1
1	line1	mRNA	1	36	.	+	0	ID=T3;Parent=G1
1	line1	gene	1	36	.	+	0	ID=G1
1	line1	exon	10	20	.	+	0	ID=E4;Parent=T4
1	line1	exon	27	32	.	+	0	ID=E5;Parent=T4
1	line1	mRNA	10	32	.	+	0	ID=T4;Parent=G2
1	line1	gene	10	32	.	+	0	ID=G2
__SAM__
read1	0	1	3	255	2S3M
read2	0	1	4	255	4M5H
read3	0	1	6	255	4M21N3M
read4	0	1	3	255	4M8N5M9N3M
read5	0	1	3	255	4M25N3M
read6	0	1	16	255	5M6N4M
