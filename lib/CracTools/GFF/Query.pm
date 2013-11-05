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

  CracTools::GFF::Query - Query GFF files easily.

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

  CracTools::GFF::Query is a tool to query GFF files without building a database.
  It is memory efficient and designed to run fast.
  You can easily retrives GFF data from a specific region of position.
  This tool can be use with CracTools::GFF::Annotation in order to parse GFF line
  into a nice usable Perl Object.

=cut

package CracTools::GFF::Query;

use strict;
use warnings;

#use Storable; # for persistency

use Set::IntervalTree;
use Fcntl qw( SEEK_SET );
use Carp;

=head1 METHODS

=head2 new

  Arg [1] : String - GFF file

  Example     : my $gffQuery = CracTools::GFF::Query->new('annotations.gff');
  Description : Create a new GFF Query object
  ReturnType  : CracTools::GFF::Query
  Exceptions  : none

=cut

sub new {
  my $class = shift;

  my $gff_file = shift;

  my $self = bless {
    GFF_FILE => $gff_file,
  }, $class;

  $self->_init();

  return $self;
}

=head2 fetchByRegion

  Arg [1] : String $seq_region_name
            The name of the sequence region that the slice will be
            created on.
  Arg [2] : int $start
            The start of the slice on the sequence region
  Arg [3] : int $end
            The end of the slice on the sequence region
  Arg [4] : int $strand
            The orientation of the slice on the sequence region

  Example     : my @annotations = $gffQuery->fetchByRegion('1',298345,309209,'+');
  Description : Retrives GFF lines that belong to the region.
  ReturnType  : Reference to an Array of strings
  Exceptions  : none

=cut

sub fetchByRegion {
  my ($self,$chr,$pos_start,$pos_end,$strand) = @_;
  my $annotations_ref = $self->_getAnnotations($chr,$strand);

  # pos_start -1 beacause Interval tree use [a,b) intervals
  my $seek_values = $self->_getAnnotations($chr,$strand)->fetch($pos_start-1,$pos_end);

  my $gff_fh = $self->_gffFilehandle;

  my @gff_lines;
  foreach (@$seek_values) {
    seek($gff_fh,$_,SEEK_SET);
    my $annot = <$gff_fh>;
    push(@gff_lines,$annot);
  }

  return \@gff_lines;
}

=head2 fetchByLocation

  Arg [1] : String $seq_region_name
            The name of the sequence region that the slice will be
            created on.
  Arg [2] : int $position
            Location to look for
  Arg [3] : int $strand
            The orientation of the slice on the sequence region

  Example     : my @annotations = $gffQuery->fetchByLocation('1',298345,'+');
  Description : Retrives GFF lines that belong to the location.
  ReturnType  : Reference to an Array of strings
  Exceptions  : none

=cut

sub fetchByLocation {
  my ($self,$chr,$position,$strand) = @_;
  return $self->fetchByRegion($chr,$position,$position,$strand);
}

sub gffFile {
  my $self = shift;
  return $self->{GFF_FILE};
}

sub _gffFilehandle {
  my $self = shift;
  return $self->{gff_fh};
}

sub _init {
  my $self = shift;
  my %annotations;

  if(!defined $self->{GFF_FILE}) {
    confess "Missing GFF file argument";
  }

  open(IN,$self->{GFF_FILE}) or die ("Cannot open file ".$self->{GFF_FILE});

  my $curr_pos = tell(IN);
  while(<IN>) {
    # skip headers
    if($_ =~ /^#/) {
      next;
    }
    my $pos = $curr_pos;
    my ($chr,$source,$feature,$start,$end,$score,$strand) = split("\t",$_,8);
    if(defined $start && defined $end && defined $chr && defined $strand) {
      $strand = convertStrand($strand);
      my $key = $self->_getAnnotationHashKey($chr,$strand);
      if(!defined $annotations{$key}) {
        $annotations{$key} = Set::IntervalTree->new;
      }
      # Minus one because gff is 1-based
      $annotations{$key}->insert($pos,$start-1,$end);
    }
    $curr_pos = tell(IN);
  }
  close IN;

  my $gff_fh;
  open($gff_fh,$self->{GFF_FILE}) or die ("Cannot open file ".$self->{GFF_FILE});

  $self->{gff_fh} = $gff_fh;
  $self->{ANNOTATIONS} = \%annotations;

}

sub _getAnnotationHashKey {
  my ($self,$chr,$strand) = @_;
  return "$chr"."@"."$strand";
}

sub _extractAnnotationHashKey {
  my $self = shift;
  my $key = shift;
  my ($chr,$strand) = split("@",$key);
  return ($chr,$strand);
}

sub _getAnnotations {
  my ($self,$chr,$strand) = @_;
  return $self->{ANNOTATIONS}{$self->_getAnnotationHashKey($chr,$strand)};
}


sub convertStrand($) {
  my $strand = shift;
  my %conversion_hash = ( '+' => 1, '-' => -1, 1 => '+', -1 => '-');
  return $conversion_hash{$strand};
}

1;
