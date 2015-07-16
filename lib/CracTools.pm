package CracTools;
# ABSTRACT: A set of tools designed to extract data from CRAC's SAM files and to provide annotations.

our $PACKAGE_NAME = "CracTools";

=head1 DESCRIPTION

CracTools-core is the cornerstone of the CracTools. It is a toolbox that aim to
ease the build of pipelines in the field of bioinformatics. It has been
originally built to produce pipelines on top of
L<CRAC|http://crac.gforge.inria.fr/> software, but you can use the
CracTools-core tools in an other context.  It has a lot of built-in features to
parse file, intersect biological events, integrate annotation, sharing
configuration.

CracTools-core is also shiped with some binaries that are directly based on the CracTools-core API:

=over 1

=item C<cractools extract>: this tools aims to extract biological events (splices, 
  snp, indels, chimeras) from BAM files produced by  CRAC's analysis.

=item C<cractools gtf2togff3>: is a tools that convert gtf annotation files to gff3
  format that is the standard in CracTools.

=item More tools are about to come (soon)

=back

=head1 MODULES

=head2 File parsing

=head3 L<CracTools::Utils>

Is a module that provide usefull functions for opening files (I/O) with iterators, 
simple parsing of standard files format (VCF,BED,GTF,GFF), or performing transormations
like reverse-complementing.

=head3 L<CracTools::SAMReader> and L<CracTools::SAMReader::SAMline>

Are modules that provide iterators and objects to easily read SAM/BAM file
generated by CRAC and provide dedicated methods to extract additional fields
added by CRAC to each record.

=head3 L<CracTools::GFF::Annotation>

Is a module to parse and access GFF3 file.

=head2 Genomic-based datastructures

=head3 L<CracTools::Interval::Query>

Is a module to store and query variables associated with genomic intervals. It
is based on the interval tree datastructure provided by L<Set::IntervalTree>.

=head3 L<CracTools::Interval::Query::File>

Acts like L<CracTools::Interval::Query> but read interval from files and return
lines of the file matching the query. It has built-in methods to parse, SAM,
G{T|F}F, BED, VCF files but you can provide your own method for other file
formats.

=head3 L<CracTools::GenomeMask>

Is a module that define a BitVector mask over a whole genome and provide method
to query this mask. It can read genome sequence and length from various sources
(SAM headers, CRAC index, User input).

=head2 Annotation

=head3 L<CracTools::Annotator>

Is a module based on L<CracTools::Interval::Query::File> that provides powerfull
methods to query annotation files and prioritize hits to fit specific
application needs.

=head2 Utilities

=head3 L<CracTools::Config>

Is a module that aim to integrate a common configuration file among all the
cractools pipelines. It automatically load the configuration file by looking to
diverse locations, then it provides methods to retrieved the variables declared
in the configuration file.

=head3 L<CracTools::Output>

Is a module that provide methods to write customized column-based output files
with pre-defined headers.

=cut

1; 
