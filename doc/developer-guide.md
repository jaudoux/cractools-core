### Accessing the options in the config file

### Opening and writing files

### Parsing a file

Parsing a file is a very easy task after you have decalre `use CracTools::Utils` on top of your scripts.

Parsing file with CracTools, is very fast since we do not built complex object, bwe only store informations in handy hash structures. See the benchmark to compare effectiveness of our solution in contrast of Bio::Perl

In order to parse files (i.e extract and format information), you can ask
CracTools-core to give you an iterator object that will help you iterate on
every elements that is present in the file. Various type of file are already supported by the CracTools-core : [S|B]AM, FAST[A|Q], VCF, BED, G[T|F]F

#### Standard files : [S|B]AM, VCF, BED, G[T|F]F

If your file is encoded in one of this format you can directly ask CracTools::Utils for your iterator

  my $it = CracTools::Utils->getFileIterator(file => myfile.vcf, type => 'vcf');

Once you have an iterator, you can iterate on all the records that are stored in the file with `while()` loop :

    while(my $vcf_line = $it->()) {
      my $ref_length = length $vcf_line->{ref}
      foreach my $alt (@{$vcf_line->{alt}}) {
        my $alt_length = length $alt;
        if($ref_length > $alt_length) {
          say "This is an deletion";
        } elseif($ref_length < $alt_length) {
          say "This is an insertion";
        } else {
          say "This is a substitution";
        }
      }
    }

#### Sequence files : FASTA, FASTQ

Sequence files can be accessed with iterators as well :

    my $seq_it = CracTools::Utils->seqFileIterator(file,"fasta"); 

You can the iterate on the "blocks" of sequence with that iterator
  
    while(my $entry = $seq_it->()) {
      print "Sequence name   : $entry->{name}
             Sequence        : $entry->{seq}
             Sequence quality: $entry->{qual}","\n";
    }

Note: if you are using paired-end sequence file, there is a special iterator that helps you iterate on the both sequence files :

    my $it = pairedEndSeqFileIterator($file);
    while (my $entry = $it->()) {
      print "Read_1 : $entry->{read1}->{seq}
             Read_2 : $entry->{read2}->{seq}";
    }


#### Other files

If you want to parse a file and the format is not yet supported, you just need to tell the CracTools how your document is formated.
In order to do that, create a function in your script that will receive a line of pure texte and return a structures hash reference.

For example we have a file with two fields


    sub parseMyFile {
      my $line = shift;
      my @fields = split(";",$line);
      return {
        chr => $fields[1],
        region => $fields[2],
        ....
      }
    }

Then you ask the CracTools for a personalized iterator, specify the argument
`parsing_method` with the subroutine you have declare in order to get the
iterator parsing you file. You can also specify a `header_regex` to skip header lines based
on a patter (for example lines beginning with a @).

    my $it = getFileIterator(file => 'myFile.csv',
      parsing_method => \$parseMyFile,
      header_regex   => '^@',
    );

And you are ready to go!

### Store information associated to genomic regions

There is a lots of case when you want to store a piece of information that is associated to a given genomic region.
Then you want to have a fast access to the information you have stored by querying a genomic region.

This *use case* can be easily achieved with CracTools::Interval::Query object.

First you have to create a new query object :

    my $interval_query = CracTools::Interval::Query->new();

Then you can add a data associated to a genomic region. You can add any type of data as long as it is a scalar.
Thus, it can be a number, a string, or an array/hash reference.

    $interval_query->addInterval($chr,$start,$end,$strand,$data);

Start position should be inferior to the end position, and strand is encoded in (1,-1) format.

After you have added some regions to the `$interval_query` object, you can query a genomic region and retrieve
the data have been added with a region that overlaps the querying region.

    my $hits = $interval_query->fetchByRegion($chr,$position,$strand)

    foreach my hit (@{$hits}) {
      say "Something...":
    }

You can also find an the closest upstream and downstream regions with `fetchNearestDown` and `fetchAllNearestUp`

### Retrieving lines in a file

In order to search lines in

### Intersect two files

### Provide annotation on a biological event

### Reading a CRAC SAM file
