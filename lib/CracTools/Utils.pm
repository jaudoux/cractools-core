###############################################################################
#                                                                             #
#    Copyright © 2012-2013 -- IRB/INSERM                                      #
#                            (Institut de Recherche en Biothérapie /          #
#                             Institut National de la Santé et de la          #
#                             Recherche Médicale)                             #
#                                                                             #
#  Auteurs/Authors:  Jerôme AUDOUX <jerome.audoux@univ-montp2.fr>             #
#                    Nicolas PHILIPPE <nicolas.philippe@inserm.fr>            #
#                                                                             #
#  -------------------------------------------------------------------------  #
#                                                                             #
#  Ce fichier fait partie de la suite CracTools qui contient plusieurs pipeline# 
#  intégrés permettant de traiter les évênements biologiques présents dans du #
#  RNA-Seq. Les CracTools travaillent à partir d'un fichier SAM de CRAC et d'un# 
#  fichier d'annotation au format GFF3.                                       #
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
#  This file is part of the CracTools which provide several integrated        #
#  pipeline to analyze biological events present in RNA-Seq data. CracTools   #
#  work on a SAM file generated by CRAC and an annotation file in GFF3 format.#
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

package CracTools::Utils;
# ABSTRACT: A set of useful functions

use strict;
use warnings;

use Carp;

=head2 reverseComplement

  Arg [1] : String - a DNA sequence
  Example : $reverse = reverseComplement('ATGC');
  Description : Reverse complemente the sequence in argument.
                For example : reverse_comptement('ATGC') returns : 'GCAT'.
  ReturnType  : String
  Exceptions  : none

=cut

sub reverseComplement($) {
  my $dna = shift;
  
  # reverse the DNA sequence
  my $revcomp = reverse $dna;
  
  # complement the reversed DNA sequence
  $revcomp =~ tr/ACGTacgt/TGCAtgca/;
  return $revcomp;
}

=head2 reverse_tab

  Arg [1] : String - a string with values separated with coma.
  Example : $reverse = reverse_tab('2,1,1,1,0,0,1');
  Description : Reverse the values of the string in argument.
                For example : reverse_tab('1,2,0,1') returns : '1,0,2,1'.
  ReturnType  : String
  Exceptions  : none

=cut

sub reverse_tab($) {
  my $string = shift;
  my @tab = split(/,/,$string);
  my $newString;
  if(@tab > 0) {
    for (my $i=$#tab ; $i > 0 ; $i--){
      $newString .= $tab[$i];
      $newString .= ",";
    }
    $newString .= $tab[0];
  }
  return $newString;
}

=head2 convertStrand
  
  Arg [1] : Character - strand

  Description : Retrun the strand converted + => 1 or 1 => +
=cut

my %conversion_hash = ( '+' => 1, '-' => '-1', '1' => '+', '-1' => '-');
sub convertStrand($) {
  my $strand = shift;
  return $conversion_hash{$strand};
}

1;
