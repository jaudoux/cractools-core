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

package CracTools::Config;

use strict;
use warnings;
use POSIX;
use utf8;

use Config::FileManager;
use Config::Simple;
use File::Basename;
use Carp;

use CracTools;

require Exporter;
our @ISA = qw(Exporter);

our $CONFIG_NAME = $CracTools::PACKAGE_NAME.".cfg";

our @EXPORT = qw(LoadConfig PrintVersion getConfVar);

our %config;

#load config
our $cfg = new Config::FileManager(
    "toolname" => $CracTools::PACKAGE_NAME, # Mandatory
    "version" => "$CracTools::VERSION", # Not mandatory
    "filename" => "$CONFIG_NAME", # Not mandatory
    "paths" => [
		".",
		".".$CracTools::PACKAGE_NAME,
		"__APPDIR__",
		"/usr/local/etc/".$CracTools::PACKAGE_NAME,
		"/etc/".$CracTools::PACKAGE_NAME
	       ], # Not mandatory
    "interactive" => 1, # Not mandatory
    );

my $default_content = "# Default configuration file __VERSION__\n#\n\n";

$cfg->defaultContent($default_content);

=head2 PrintVersion

Print (in an uniformized way) the version information of the CracUtil script.

Usage:

  PrintVersion();

=cut

sub PrintVersion() {
  printf( "Script '%s' from %s v. %s (%s v. %s)\n",
	  basename($0),
	  $CracTools::PACKAGE_NAME, $CracTools::VERSION,
	  $CracTools::PACKAGE_NAME, $CracTools::VERSION);
}

=head2 LoadConfig

Usage:

  my $default_cfg = LoadConfig(); # Use (and get) default config file (see L<Config::FileManager>)
  LoadConfig("your_file.conf");

=cut

sub LoadConfig(;$) {
    my ($config_file) = @_;
    if (!defined $config_file) {
	$cfg->update();
	$config_file = $cfg->getPath();
    }
    Config::Simple->import_from($config_file, \%config);
    return $config_file;
}

=head2 getConfVar

Usage:

  my $var = getConfVar("GENOME"); 

=cut

sub getConfVar(;$) {
  my $var_name = shift;
  my $die = shift;
  if(defined $config{$var_name}) {
    return $config{$var_name};
  } else {
    if(defined $die && $die eq 1) {
      croak("Config variable \"$var_name\" not found.");
    } else {
      return undef;
    }
  }
}

1; 
__END__

=head1 AUTHORS

Nicolas PHILIPPE E<lt>L<nicolas.philippe@inserm.fr|mailto:nicolas.philippe@inserm.fr>E<gt>.
Alban MANCHERON E<lt>L<alban.mancheron@lirmm.fr|mailto:alban.mancheron@lirmm.fr>E<gt>,

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
