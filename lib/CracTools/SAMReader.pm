package CracTools::SAMReader;
# ABSTRACT: An easy to use tool to read files in SAM format.

=head1 SYNOPSIS

  use CracTools::SAMReader;

  # Creating the reader
  my $sam_reader = CracTools::SAMreader->new($sam,'CRAC');

  # Get an iterator to go through the SAM file in a linear way
  my $it = $sam_reader->iterator();
  
  # Iterate on lines and explore CRAC special fields of SAM
  while(my $line = $it->()) {
    if(defined $line->events('Junction') && $line->isClassified('normal')) {
      my @junctions = @{$line->events('Junction')};
      foreach my $junction (@junctions) {
        print "Foud Junction : [type : $junction->{type}, loc : $junction->{loc}, gap : $junction->{gap}]\n";
      } 
    }
  }

=head1 DESCRIPTION

Reader for SAM format, including CRAC special fields.

=cut

use strict;
use warnings;
use Carp;
use CracTools::SAMReader::SAMline;

=head1 METHODS

=head2 new

  Arg [1] : String - SAM file
  Arg [2] : (Optional) String - SAM type
            - CRAC
            - CRAC_EMT

  Example     : $reader = CracTools::SAMreader->new('file.sam','CRAC');
  Description : Create a new reader object
  ReturnType  : CracTools::SAMreader
  Exceptions  : none

=cut

sub new {
  my $class = shift;
  my ($sam_file,$sam_type) = @_;

  my $self = bless{ 
      sam_file => $sam_file,
      sam_type => $sam_type,
  }, $class;

  $self->init();

  return $self;
}

=head2 iterator

  Example     : my $it = $sam_reader->iterator();
                while(my $line = $it->()) {
                  print $line->seq,"\n";
                }
  Description : Create an iterator to go throud each lines of the file
  ReturnType  : Iterator on CracTools::SAMline
  Exceptions  : none

=cut

sub iterator {
  my $self = shift;
  my $f_it = $self->iteratorFile("IGNORE_HEADERS");

  return sub {
    my ($line) = $f_it->();
    my $sam_line;
    if(defined $line) {
      $sam_line = CracTools::SAMReader::SAMline->new($line);
    }
    return $sam_line;
  };
}

=head2 iteratorFile

  Arg [1] : (Optional) String - options (IGNORE_HEADERS,..)

  Example     : my $it_f = $sam_reader->iteratorFile();
                while(my ($line,$line_number) = $it->()) {
                  print $line,"\n";
                }
  Description : Create an iterator to go throud each lines of the file
  ReturnType  : Iterator on Array (String,Int) where the <String> is the
                line, and <Int> the line number.
  Exceptions  : none

=cut

sub iteratorFile {
  my $self = shift;
  my $option = shift;
  my $sam_file = $self->{sam_file};

  if($sam_file =~ /\.sam$/) {
    open(SAM,"< $sam_file") or die ("Cannot open $sam_file");
  } elsif($self->{sam_file} =~ /\.sam.gz$/) {
    open(SAM,"gunzip -c $sam_file |") or die ("Cannot open $sam_file");
  } elsif($self->{sam_file} =~ /\.bam$/) {
    open(SAM, "-|", "samtools view -h $sam_file" )or die "Cannot open $sam_file, check if samtools are installed.";
  } else {
    open(SAM,"< $sam_file") or die ("Cannot open $sam_file");
    warn "Unknown file format. We assume this is SAM (uncompressed).";
  }

  my $next_line;
  my $line_number = 0;

  if(defined $option && $option eq "IGNORE_HEADERS") {
    while(my $line = <SAM>) {
      if(!($line =~ /^@/)) {
        $next_line = $line;
        $line_number++;
        last;
      }
    }
  } else {
    $next_line = <SAM>;
  }

  return sub {
    my $sam_line = $next_line;
    $next_line = <SAM>;
    $line_number++;
    if($sam_line) {
      return $sam_line, $line_number;
    } else {
      close(SAM) or warn $! ? "Error closing samtools pipe: $!" : "Exit status $? from samtools";
      return ();
    }
  };
}

=head1 GETTERS AND SETTERS

=cut

=head2 header

  Description : Getter/setter for attribute header
  ReturnType  : none
  Exceptions  : none

=cut

sub header {
  my $self = shift;
  return $self->{header};
}

=head2 refSeqLength

  Description : Return the length of the reference sequence given in argument
  ReturnType  : Integer

=cut

sub refSeqLength {
  my $self = shift;
  my $ref_seq = shift;
  croak("Missing reference sequence name in arguement") unless defined $ref_seq;
  my $refseq_lengths = $self->allRefSeqLengths();
  return $refseq_lengths->{$ref_seq};
}

=head2 allRefSeqLengths

Return a hashref with all reference sequence names and length.

  my $refseq_lengths = $sam_reader->allRefSeqLength();

The return value is :

  { refseq_name => length,
    refseq_name => lenght,
    ...
  }

=cut

sub allRefSeqLengths {
  my $self = shift;
  my @header_lines = split('\n',$self->header);
  my %refseq_lengths; 
  foreach (@header_lines) {
    if ($_ =~/\@SQ.*SN:/) {
      my ($name,$length) = $_ =~/\@SQ.*SN:(\S+)\s+LN:(\d+)+/;    
      $refseq_lengths{$name} = $length;
    }
  }
  return \%refseq_lengths;
}

=head2 commandLine

  Description : Return crac command line defined in SAM's header
  ReturnType  : String

=cut

sub commandLine {
  my $self = shift;
  if(defined $self->header) {
    my @header_lines = split('\n',$self->header);
    my $command_line; 
    foreach (@header_lines) {
      if ($_ =~/\@PG.*PN:crac/) {
        ($command_line) = $_ =~ /CL:([^\t]+)/;    
      }
    }
    return $command_line;
  } else {
    return undef;
  }
}

=head2 getCracArgumentValue

  Description : Retrun the value of the specified argument in crac command line

=cut 

sub getCracArgumentValue {
  my $self = shift;
  my $argument = shift;
  my $command_line = $self->commandLine;
  if(defined $command_line) {
    my ($value) = $command_line =~ /--$argument\s+(\S+)/;
    return $value;
  } else {
    return undef;
  }
}

=head2 hasCracOption

  Description : Return true if crac command line has specified a certain option

=cut

sub hasCracOption {
  my $self = shift;
  my $option = shift;
  croak("Missing argument") unless defined $option;
  if(defined $self->commandLine) {
    return $self->commandLine =~ /--$option/;
  } else {
    return 0;
  }
}

=head2 getCracVersionNumber

  Description : Return CRAC version number

=cut

sub getCracVersionNumber {
  my $self = shift;
  if(defined $self->header) {
    my @header_lines = split('\n',$self->header);
    my $version_number; 
    foreach (@header_lines) {
      if ($_ =~/\@PG.*PN:crac/) {
        ($version_number) = $_ =~ /VN:([^\t]+)/;    
        last;
      }
    }
    return $version_number;
  } else {
    return undef;
  }
}

=head1 PRIVATE METHODS

=head2 init (private)

  Description : Initialization method
  ReturnType  : none
  Exceptions  : none

=cut

sub init {
  my $self = shift;
  my $f_it = $self->iteratorFile;
  my $header;
  while(my ($line) = $f_it->()) {
    if($line =~ /^@/) {
      $header .= $line;
    } else {
      last;
    }
  }
  $self->{header} = $header;
}

1;
