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

  CracTools::GFF::Annotation - Parse GFF lines.

=head1 SYNOPSIS

Usage:

  use CracTools::GFF::Query;

  # Creating the reader
  my $gffQuery = CracTools::GFF::Query->new('annotations.gff');

  my @annotations = $gffQuery->fetchByLocation('1',298345,'+');

  foreach my $gff_line (@annotations) {
    my $annotation = CracTools::GFF::Annotation->new($gff_line);
    print "Gene_id : ",$annotation->getAttribute('gene_id'),"\n";
  }

=head1 DESCRIPTION

  CracTools::GFF::Annotataion is an object to parse and access GFF line's fields.

=cut

package CracTools::GFF::Annotation;

use Carp;

=head1 METHODS

=head2 new

  Arg [1] : String - $line
            GFF line
  Arg [2] : String - $format (optional) - default 'gff2'
            GFF format (gff2 or gff3)

  Example     : my $annotation = CracTools::GFF::Annotation->new($gff_line);
  Description : Create a new CracTools::GFF::Annotation object
                If a gff line is passed in argument, the line will be parsed
                and loaded.
  ReturnType  : CracTools::GFF::Query
  Exceptions  : none

=cut

sub new {
  my $class = shift;
  my $line = shift;
  my $format = shift;
  if(!defined $format) {
    $format = 'gff2';
  }

  my $self = bless {format => $format}, $class;

  if(defined $line) {
    $self->_init($line);  
  }

  return $self;
}

sub _init {
  my ($self,$line) = @_;

  # Split the line with TABS
  my ($chr,$source,$feature,$start,$end,$score,$strand,$phase,$attributes) = split("\t",$line);

  # Get the strand if 1/-1 format
  $strand = convertStrand($strand);

  # Loading the 8 first columns
  $self->{chr} = $chr;
  $self->{source} = $source;
  $self->{feature} = $feature;

  # Conversion to 0-based coordinate system
  $self->{start} = $start-1;
  $self->{end} = $end-1;

  $self->{score} = $score;
  $self->{strand} = $strand;
  $self->{phase} = $phase;

  # Loading attributes
  my @attributes_tab = split(";",$attributes);
  foreach my $attr (@attributes_tab) {
    my ($k,$v);
    if($self->{format} =~ /gff3/i) {
      ($k,$v) = $attr =~ /(\S+)=(.*)/;
    } else {
      ($k,$v) = $attr =~ /(\S+)\s+"(.*)"/;
    }
    if(defined $k && defined $v) {
      $self->{attributes}{$k} = $v;
    } #else {
    #  carp("Error parsing attribute $attr");
    #}
  }

}

sub chr {
  my $self = shift;
  my $chr = shift;
  if(defined $chr) {
    $self->{chr} = $chr;
  }
  return $self->{chr};
}

sub source {
  my $self = shift;
  my $source = shift;
  if(defined $source) {
    $self->{source} = $source;
  }
  return $self->{source};
}

sub feature {
  my $self = shift;
  my $feature = shift;
  if(defined $feature) {
    $self->{feature} = $feature;
  }
  return $self->{feature};
}

sub start {
  my $self = shift;
  my $start = shift;
  if(defined $start) {
    $self->{start} = $start;
  }
  return $self->{start};
}

sub end {
  my $self = shift;
  my $end = shift;
  if(defined $end) {
    $self->{end} = $end;
  }
  return $self->{end};
}

sub score {
  my $self = shift;
  my $score = shift;
  if(defined $score) {
    $self->{score} = $score;
  }
  return $self->{score};
}

sub strand {
  my $self = shift;
  my $strand = shift;
  if(defined $strand) {
    $self->{strand} = $strand;
  }
  return $self->{strand};
}

sub gffStrand {
  my $self = shift;
  return convertStrand($self->{strand});
}

sub phase {
  my $self = shift;
  my $phase = shift;
  if(defined $phase) {
    $self->{phase} = $phase;
  }
  return $self->{phase};
}

sub parents {
  my $self = shift;
  if(defined $self->attribute('Parent')) {
    return split(',',$self->attribute('Parent'));
  } else {
    return ();
  }
}

sub attribute {
  my $self = shift;
  my $key = shift;
  my $value = shift;
  if(defined $value) {
    return $self->{attributes}{$key} = $value;
  } elsif(defined $key) {
    return $self->{attributes}{$key};
  } else {
    return undef;
    #croak ("Missing attribute key to retreive attribute value");
  }
}

sub convertStrand($) {
  my $strand = shift;
  my %conversion_hash = ( '+' => 1, '-' => -1, 1 => '+', -1 => '-');
  return $conversion_hash{$strand};
}

1;
