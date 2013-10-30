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
#                   Alban MANCHERON  <alban.mancheron@lirmm.fr>               #
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

CracTools::Output - A module to manage output files.

=head1 SYNOPSIS

  # Creating a default output object.
  # Everything will be print to the standard output
  my $output = CracTools::Output->new();
  
  # Print nice headers
  my $output->printHeaders(version => '1.01',summary => 'blabla', args => @ARGV);
  
  # This will print "foo\tbar\n"
  $output->printLine('foo','bar');
  
  # Using semicolon as separator charcater
  my $output = CracTools::Output->new(sep => ';');
  
  # Print into a file
  my $output = CracTools::Output->new(file => 'foo.bar');

=head1 DESCRIPTION

CracTools::Output is a simple tool to generate Char-Separated files.

=cut

package CracTools::Output;

use strict;
use warnings;

use File::Basename;
use CracTools;

use constant DEFAULT_SEP => "\t";

use strict;
use warnings;

sub new {
  my $class = shift;

  my %args = @_;

  my $sep = DEFAULT_SEP unless $args{sep};
  my $output = $args{file};
  my $out_stream;

  if(defined $output) {
    open($out_stream,">$output") or die ("Enable to open $output file.\n");
  } else {
    $out_stream = \*STDOUT;
  }
  
  my $self = bless {
    sep => $sep, 
    out_stream => $out_stream,
  }, $class;
  
  return $self
}

sub printHeaders {
  my $self = shift;
  my %args = @_;

  my $version = $args{version};
  my $summary = $args{summary};
  my @arguments;
  if(defined $args{args}) {
    @arguments = @{$args{args}};
  }

  $self->printlnOutput("# Date: ".localtime);
  $self->printlnOutput("# Module: $CracTools::PACKAGE_NAME (v $CracTools::VERSION)");
  if(defined $version) {
    $self->printlnOutput("# Script: ".basename($0)." (v $version)");
  } else {
    $self->printlnOutput("# Script: ".basename($0));
  }
  if(@arguments > 0) {
    $self->printOutput("# Args: ");
    foreach my $arg (@arguments) {
      $self->printOutput(" $arg");
    }
    $self->printlnOutput();
  }
  if(defined $summary) {
    $self->printlnOutput("# Summary:");
    my @lines = split /\n/, $summary;
    foreach my $line (@lines) {
      $self->printlnOutput("# $line");
    }
  }
  #$self->printlnOutput();
}

sub printHeaderLine {
  my $self = shift;
  my $first_field = shift;
  $first_field = '' unless defined $first_field;
  $self->printLine("# $first_field",@_);
}

sub printLine {
  my $self = shift;
  for(my $cpt = 0; $cpt < scalar @_; $cpt++) {
    if(!defined $_[$cpt]) {
      $_[$cpt] = 'NONE';
    }
  }
  $self->printlnOutput(join($self->{sep},@_));
}

sub printOutput {
  my $self = shift;
  my $stuff = shift;
  my $stream = $self->{out_stream};
  print $stream $stuff;
}

sub printlnOutput {
  my $self = shift;
  my $stuff = shift;
  if(defined $stuff) {
    $self->printOutput("$stuff\n");
  } else {
    $self->printOutput("\n");
  }
}

1;
