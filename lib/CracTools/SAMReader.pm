###############################################################################
#                                                                             #
#    Copyright © 2012-2013 -- IRB/INSERM                                      #
#                            (Institut de Recherche en Biothérapie /          #
#                             Institut National de la Santé et de la          #
#                             Recherche Médicale)                             #
#                             LIRMM/UM2                                       #
#                            (Laboratoire d'Informatique, de Robotique et de  #
#                             Microélectronique de Montpellier /              #
#                             Université de Montpellier 2)                    #
#                                                                             #
#  Auteurs/Authors: Nicolas PHILIPPE <nicolas.philippe@inserm.fr>             #
#                   Jerome AUDOUX  <jerome.audoux@univ-montp2.fr>             #
#                                                                             #
#  -------------------------------------------------------------------------  #
#                                                                             #
#  Ce fichier  fait partie  du Pipeline  de traitement  de données NGS de la  #
#  plateforme ATGC labélisée par le GiS IBiSA.                                #
#                                                                             #
#  Ce logiciel est régi  par la licence CeCILL  soumise au droit français et  #
#  respectant les principes  de diffusion des logiciels libres.  Vous pouvez  #
#  utiliser, modifier et/ou redistribuer ce programme sous les conditions de  #
#  la licence CeCILL  telle que diffusée par le CEA,  le CNRS et l'INRIA sur  #
#  le site "http://www.cecill.info".                                          #
#                                                                             #
#  En contrepartie de l'accessibilité au code source et des droits de copie,  #
#  de modification et de redistribution accordés par cette licence, il n'est  #
#  offert aux utilisateurs qu'une garantie limitée.  Pour les mêmes raisons,  #
#  seule une responsabilité  restreinte pèse  sur l'auteur du programme,  le  #
#  titulaire des droits patrimoniaux et les concédants successifs.            #
#                                                                             #
#  À  cet égard  l'attention de  l'utilisateur est  attirée sur  les risques  #
#  associés  au chargement,  à  l'utilisation,  à  la modification  et/ou au  #
#  développement  et à la reproduction du  logiciel par  l'utilisateur étant  #
#  donné  sa spécificité  de logiciel libre,  qui peut le rendre  complexe à  #
#  manipuler et qui le réserve donc à des développeurs et des professionnels  #
#  avertis  possédant  des  connaissances  informatiques  approfondies.  Les  #
#  utilisateurs  sont donc  invités  à  charger  et  tester  l'adéquation du  #
#  logiciel  à leurs besoins  dans des conditions  permettant  d'assurer  la  #
#  sécurité de leurs systêmes et ou de leurs données et,  plus généralement,  #
#  à l'utiliser et l'exploiter dans les mêmes conditions de sécurité.         #
#                                                                             #
#  Le fait  que vous puissiez accéder  à cet en-tête signifie  que vous avez  #
#  pris connaissance  de la licence CeCILL,  et que vous en avez accepté les  #
#  termes.                                                                    #
#                                                                             #
#  -------------------------------------------------------------------------  #
#                                                                             #
#  This File is part of the NGS data processing Pipeline of the ATGC          #
#  accredited by the IBiSA GiS.                                               #
#                                                                             #
#  This software is governed by the CeCILL license under French law and       #
#  abiding by the rules of distribution of free software. You can use,        #
#  modify and/ or redistribute the software under the terms of the CeCILL     #
#  license as circulated by CEA, CNRS and INRIA at the following URL          #
#  "http://www.cecill.info".                                                  #
#                                                                             #
#  As a counterpart to the access to the source code and rights to copy,      #
#  modify and redistribute granted by the license, users are provided only    #
#  with a limited warranty and the software's author, the holder of the       #
#  economic rights, and the successive licensors have only limited            #
#  liability.                                                                 #
#                                                                             #
#  In this respect, the user's attention is drawn to the risks associated     #
#  with loading, using, modifying and/or developing or reproducing the        #
#  software by the user in light of its specific status of free software,     #
#  that may mean that it is complicated to manipulate, and that also          #
#  therefore means that it is reserved for developers and experienced         #
#  professionals having in-depth computer knowledge. Users are therefore      #
#  encouraged to load and test the software's suitability as regards their    #
#  requirements in conditions enabling the security of their systems and/or   #
#  data to be ensured and, more generally, to use and operate it in the same  #
#  conditions as regards security.                                            #
#                                                                             #
#  The fact that you are presently reading this means that you have had       #
#  knowledge of the CeCILL license and that you accept its terms.             #
#                                                                             #
###############################################################################

=head1 NAME

  CracTools::SAMReader - An easy to use tool to read files in SAM format.

=head1 SYNOPSIS

Usage:

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

package CracTools::SAMReader;

use strict;
use warnings;
use CracTools::SAMReader::SAMline;

=head1 METHODS

=head2 new

  Arg [1] : String - SAM file
  Arg [2] : (Optional) String - SAM type
            - CRAC
            - CRAC_EMT

  Example     : $reader = CracTools::SAMreader->new('file.sam','CRAC');
  Description : Create a new reader obect
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
  my $f_it = $self->iteratorFile("INGNORE_HEADERS");

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

  Arg [1] : (Optional) String - options (INGORE_HEADERS,..)

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
    open(SAM, "-|", "samtools view $sam_file" )or die "Cannot open $sam_file, check if samtools are installed.";
  } else {
    die "Unknown file format. Must be either a BAM or a SAM(.gz)";
  }

  my $next_line;
  my $line_number = 0;

  if(defined $option && $option eq "INGNORE_HEADERS") {
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

sub commandLine {
  my $self = shift;
  my @header_lines = split('\n',$self->header);
  my $command_line; 
  foreach (@header_lines) {
    if ($_ =~/\@PG.*PN:crac/) {
      ($command_line) = $_ =~ /CL:([^\t]+)/;    
    }
  }
  return $command_line;
}

# retrun the value of the specified argument in crac command line
sub getCracArgumentValue {
  my $self = shift;
  my $argument = shift;
  my $command_line = $self->commandLine;
  my ($value) = $command_line =~ /--$argument\s+(\S+)/;
  return $value;
}

# returne true if crac command line has specified a certain option
sub hasCracOption {
  my $self = shift;
  my $option = shift;
  croak("Missing argument") unless defined $option;
  return $self->commandLine =~ /--$option/;
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

=head1 AUTHORS

Jerome AUDOUX E<lt>L<jerome.audoux@etud.univ-montp2.fr|mailto:jerome.audoux@univ-montp2.fr>E<gt>.

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2012-2013 -- IRB/INSERM
                           (Institut de Recherche en Biothérapie /
                            Institut National de la Santé et de la
                            Recherche Médicale)
                           LIRMM/UM2
                           (Laboratoire d'Informatique, de Robotique et de
                            Microélectronique de Montpellier /
                            Université de Montpellier 2)

=head2 FRENCH

Ce fichier  fait partie  du Pipeline  de traitement  de données NGS de la
plateforme ATGC labélisée par le GiS IBiSA.

Ce logiciel est régi  par la licence CeCILL  soumise au droit français et
respectant les principes  de diffusion des logiciels libres.  Vous pouvez
utiliser, modifier et/ou redistribuer ce programme sous les conditions de
la licence CeCILL  telle que diffusée par le CEA,  le CNRS et l'INRIA sur
le site "http://www.cecill.info".

=head2 ENGLISH

This File is part of the NGS data processing Pipeline of the ATGC
accredited by the IBiSA GiS.

This software is governed by the CeCILL license under French law and
abiding by the rules of distribution of free software. You can use,
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info".

=cut

1;
