=head1

    Crac Wrapper

=cut

package CracTools::Aligner::Crac;
use CracTools::SAMReader;
use strict;
use warnings;
use File::Temp;
use Carp;
use Data::Dumper;
use CracTools::Const;

sub new {

  my $class = shift;
  my @args = @_;
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
  my $command = " ";  
  #my $command=print MAPPER," ", print KMER," ",print INDEX;
  my $sam_file = new File::Temp( SUFFIX => '.sam');
  if(@args) {  
    $command=$command." -r ";
    foreach(@args) {
        $command=$command." ".$_; 
    }
  }
  $command=$command." -o ".$sam_file;
  #Run CRAC
  system($command) == 0 or die "Can't execute $command";
   
  #Output
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
