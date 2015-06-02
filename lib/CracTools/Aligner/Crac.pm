=head1

    Crac Wrapper

=cut

package CracTools::Aligner::Crac;
use CracTools::SAMReader;
use CracTools::SAMReader::SAMline;
use strict;
use warnings;
use File::Temp;
use Carp;
use Data::Dumper;

#Default
our $mapper = "/data/projects/crac-dev/src/crac";
my $index = "-i /data/indexes/crac/GRCh38";
my $k = "-k 20";
my $options = "--detailed-sam";

sub new {

  my $class = shift;
  my @args = @_;
  my $command=$mapper;
  my $self = bless {
    samList => []
  }, $class;
  if (@args) {
    $self->_init(@args); 
  } 
  return $self;
}

sub _init {

  my $self = shift;
  my @args = @_;
  my $command = $mapper;
  my $nb_alignements;
  
  #Command line
  foreach(@args) {
    $command = $command." ".$_;
  }
  my $sam_file = new File::Temp( SUFFIX => '.sam');
  $command=$command." -o ".$sam_file;
  
  #Run CRAC
  system($command) == 0 or die "Can't execute $command";
  
  #Output parsing
  my $sam_reader = CracTools::SAMReader->new($sam_file);
  close $sam_file;
  my $it = $sam_reader->iterator();
  my @list_align;
  while(my $line = $it->()) {
    push(@{$self->{samList}},$line);
  }
  return $self;  
}

1;
