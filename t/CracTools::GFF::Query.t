#! /usr/bin/perl
#
use strict;
use warnings;

use Test::More tests => 16;
use CracTools::GFF::Query;
use File::Temp;

my $line1 = "1\tEnsembl\texon\t1\t10\t.\t+\t0\n";
my $line2 = "1\tEnsembl\texon\t2\t3\t.\t+\t0\n";
my $line3 = "1\tEnsembl\texon\t3\t9\t.\t+\t0\n";
my $line4 = "1\tEnsembl\texon\t1\t2\t.\t-\t0\n";
my $line5 = "2\tEnsembl\texon\t4\t8\t.\t+\t0\n";

my $gff = $line1.$line2.$line3.$line4.$line5;


# Create a temp file with the SAM lines described above
my $gff_file = new File::Temp( SUFFIX => '.sam', UNLINK => 1);
print $gff_file $gff;
close $gff_file;

my $gffQuery = CracTools::GFF::Query->new($gff_file);

# Testing private methods
ok($gffQuery->_extractAnnotationHashKey($gffQuery->_getAnnotationHashKey(1,'1')) ~~ (1,'1'), 'Extract annotations keys');

ok($gffQuery->_getAnnotationHashKey(1,'1') ne $gffQuery->_getAnnotationHashKey(1,'-1'), 'Unicity of annotations keys (1)');
ok($gffQuery->_getAnnotationHashKey(1,'1') ne $gffQuery->_getAnnotationHashKey(2,'1'), 'Unicity of annotations keys (2)');

ok(@{$gffQuery->_getAnnotationKeysSorted(1,'1','start')} ~~ [0,1,2], 'getAnnotationKeysSorted()');
ok(@{$gffQuery->_getAnnotationKeysSorted(1,'1','end')} ~~ [1,2,0], 'getAnnotationKeysSorted()');

ok(@{$gffQuery->_getAnnotationKeysSubset(1,'1','start',2,'lower')} ~~ [0,1], 'getAnnotationKeysSubset() (1)');
ok(@{$gffQuery->_getAnnotationKeysSubset(1,'1','start',3,'upper')} ~~ [2], 'getAnnotationKeysSubset() (2)');
ok(@{$gffQuery->_getAnnotationKeysSubset(1,'1','end',3,'lower')} ~~ [1], 'getAnnotationKeysSubset() (3)');
ok(@{$gffQuery->_getAnnotationKeysSubset(1,'1','end',4,'upper')} ~~ [2,0], 'getAnnotationKeysSubset() (4)');
ok(@{$gffQuery->_getAnnotationKeysSubset(1,'-1','end',3,'upper')} ~~ [], 'getAnnotationKeysSubset() (5)');
ok(@{$gffQuery->_getAnnotationKeysSubset(2,'1','start',3,'lower')} ~~ [], 'getAnnotationKeysSubset() (6)');

ok(@{$gffQuery->fetchByLocation(1,3,'1')} == 3, 'findAnnotations() (1)');
ok(@{$gffQuery->fetchByLocation(1,3,'-1')} == 0, 'findAnnotations() (2)');
ok(@{$gffQuery->fetchByLocation(1,4,'1')} == 2, 'findAnnotations() (3)');
ok(@{$gffQuery->fetchByLocation(1,10,'1')}[0] eq $line1, 'findAnnotations() (4)');
is(@{$gffQuery->fetchByRegion(1,3,3,1)}[0],@{$gffQuery->fetchByLocation(1,3,1)}[0],'fetchByRegion');

