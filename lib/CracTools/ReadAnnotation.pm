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

ReadAnnotation - A tool for easy annotation of a reads on the Ensembl API.

=head1 SYNOPSIS

    annoteTagsOnEnsembl.pl -g genome -i mapping_file -k gene_point -o output -l factor_length -m max_query_simultaneous [--config file.cfg] [--registry PATH_TO_EnsemblRegistry.pm] [--verbose] |--version]
    -g genome : genome which is associated with tags.
             'Homo Sapiens' for Human (default),
             'Pan troglodytes' for Chimpanzee,
             'Mus musculus' for Mouse,
             'Macaca mulatta' for Macaque,
             'Pongo pygmaeus' for Orangutan,
              etc (cf http://www.ensembl.org/info/about/species.html)
    -i mapping_file : a file with mapping information for each tag, 
                      one tag by line and one location by tag (format tag|chr,strand  position).    
    -l factor_length : the tag length. 
    -k gene_point : the distance to find the first gene in the 5' neighborhood and the 3' neighborhood of the tag (DISTANCE_MAX by default).
    -o output : output file.
    -m max_query_simultaneous : number of simultaneous queries to run on Ensembl GB database

    --config file.cfg : use file.cfg instead of default config file.
    --registry PATH_TO_EnsemblRegistry.pm
    --verbose (optional) : we can use this option to print more details of the procedure in the log file.
    --version : Ensembl API version

=head1 OPTIONS

    output :
              output                  : a csv file with tabulated format for the annotation features for each tag.

           Each line contains two level of information: 
             
              1) In case of a tag cross a protein_coding gene or a pseudogene, a relative annotation by tag is saved: 
                 3PRIM_UTR_sense      : it is in a 3PRIM_UTR of a transcript.  
                 3PRIM_UTR_antisense  : it is antisense of a 3PRIM_UTR of a transcript.  
                 CDS_sense            : it is in a CDS of a transcript.  
                 CDS_antisense        : it is antisense of a CDS of a transcript.  
                 5PRIM_UTR_sense      : it is in a 5PRIM_UTR of a transcript.  
                 5PRIM_UTR_antisense  : it is antisense of a 5PRIM_UTR of a transcript.  
                 INXON_sense          : it is overlap an exon and an intron.
                 INXON_antisense      : it is antisense of an overlaping between an exon and an intron.
                 INTRON_sense         : it is in an INTRON of a transcript.    
                 INTRON_antisense     : it is antisense of an INTRON of a transcript.  
                 INTER_PROXIMAL       : it is in an intergenic region but near of a 5' gene. 
                                        (before the gene_point value) and it can be considered as a 3'UTR variant.  
                 INTER_DISTAL         : it is in an intergenic region and it can be considered as a new transcript.  
                 INTER_DISTAL_EST     : it is intergenic distal but some EST overlap the tag.
    
              2) In case of a tag cross a non_coding gene, a relative annotation by tag is saved:
                 small_ncRNA          : it is included inside a small non_coding RNA (miRNA, snRNA, snoRNA, rRNA, Mt_rRNA, Mt_tRNA
                                                                                 , tRNA, scRNA, ncRNA, 3prime_overlapping_ncrna)
                 lincRNA              : it is included inside a long intergenic non_coding RNA
                 other_lncRNA         : it is included inside an other long non_coding RNA (antisense, sense_intronic, processed_transcript)
                 other_noncodingRNA   : it is included inside an other non_coding RNA (misc_RNA, ncRNA_host, sense_overlapping
                                                                                      , retained_intron, processed_pseudogene
                                                                                      , unprocessed_pseudogene
                                                                                      , transcribed_processed_pseudogene
                                                                                      , transcribed_unprocessed_pseudogene
                                                                                      , retrotransposed, unitary_pseudogene)


=head1 DESCRIPTION

An easy tool to annotate reads using the Ensembl API.

=cut

package CracTools::ReadAnnotation;

use strict;
use warnings;

use CracTools::Config;

use Bio::EnsEMBL::Registry;

our @type_prot = ('3PRIM_UTR_sense','CDS_sense','5PRIM_UTR_sense','INXON_sense','INTRON_sense','3PRIM_UTR_antisense','CDS_antisense','5PRIM_UTR_antisense','INXON_antisense','INTRON_antisense','INTER_PROXIMAL','INTER_DISTAL_EST','INTER_DISTAL');

our @type_noncoding = ("small_ncRNA","lincRNA","other_lncRNA","other_noncodingRNA");

our %priority_annotation = (1   =>  $type_prot[0],
                            2   =>  $type_prot[1],
                            3   =>  $type_prot[2],
                            4   =>  $type_prot[3],
                            5   =>  $type_prot[4],
                            6   =>  $type_prot[5],
                            7   =>  $type_prot[6],
                            8   =>  $type_prot[7],
                            9   =>  $type_prot[8],
                            10  =>  $type_prot[9],
                            11  =>  $type_prot[10],
                            12  =>  $type_prot[11],
                            13  =>  $type_prot[12]);


=head1 PUBLIC METHODS

=cut


=head2 new

  Arg [TAG_LENGTH]:   Integer - Tag_length to annotate (same as k in CRAC).
  Arg [SPECIES]:      (Optional) String - Genome to use for querying Ensembl
                      Default : 'Human'
  Arg [DISTANCE_MAX]: (Optional) String - Distance max constatn.
  Arg [INTERGENIC_THRESHOLD]:   (Optional) String - Gene point constant/

  Exemple     : $annotator = ReadAnnotation->new();
  Description : Create a new ReadAnnotation object
  ReturnType  : ReadAnnotation
  Exceptions  : none

=cut

sub new {
  my $class = shift;

  my ($args) = @_;

  my $registry = 'Bio::EnsEMBL::Registry';

  if(!defined $args->{TAG_LENGTH}) {
    die("Argument TAG_LENGTH not specified in ReadAnnotation.pm constructor.");
  }
 
  #my $config_file = LoadConfig($config_file);
  LoadConfig();

  my $genome =  $args->{SPECIES};
  my $distance_max =  $args->{DISTANCE_MAX};
  my $gene_point =  $args->{INTERGENIC_THRESHOLD};
  my $ensembl_registry_cfg =  $args->{ENSEMBL_REGISTRY_CFG};
  $genome = getConfVar("GENOME") unless defined $args->{SPECIES};
  $distance_max = getConfVar("DISTANCE_MAX") unless defined $args->{DISTANCE_MAX};
  $gene_point = getConfVar("INTERGENIC_THRESHOLD") unless defined $args->{INTERGENIC_THRESHOLD};
  $ensembl_registry_cfg = getConfVar("ENSEMBL_REGISTRY_CFG") unless defined $args->{ENSEMBL_REGISTRY_CFG};

  if(defined $ensembl_registry_cfg) {
    $registry->load_all($ensembl_registry_cfg) or die("Not able to load Ensembl registry conf file : $ensembl_registry_cfg\n");
  } elsif(defined $ENV{'ENSEMBL_REGISTRY'}) {
    $registry->load_all() or die("Not able to load ensembl registry from env variable ENSEMBL_REGISTRY: $ENV{ENSEMBL_REGISTRY}\n");
  } else {
    print STDERR "No Ensembl registry configuration file have been found. Connecting to Ensembl with default values.\n";
    $registry->load_registry_from_db(
        -host    => 'ensembldb.ensembl.org',
        -user    => 'anonymous',
        -verbose => '1'
    );
  }

  my $slice_adaptor = $registry->get_adaptor( $genome, 'Core', 'Slice' );
  my $sliceEST_adaptor = $registry->get_adaptor( $genome, 'OtherFeatures', 'Slice');

  my $self = bless { 
     registry => $registry,
     slice_adaptor => $slice_adaptor,
     sliceEST_adaptor => $sliceEST_adaptor,
     tag_length => $args->{TAG_LENGTH},
     species => $genome,
     DISTANCE_MAX => $distance_max,
     INTERGENIC_THRESHOLD => $gene_point,
  }, $class;

  return $self;
}



=head2 getAnnotation

  Arg [1]     : String - Chromosome
  Arg [2]     : Integer - Strand
  Arg [3]     : Integer - position
  Arg [4]     : (Optional) String ['before | 'after'] - Sense (are we looking for a tag before this position or after)
                Default : 'after'
  Arg [5]     : (Optional) String - Tag (sequence of the tag to process unit tests) 

  Exemple     : my $slice = readAnnotation->getAnnotation();
  Description : Create an annotation hash for the given position.
  ReturnType  : Annotation hash of a tag : 
                %annotation = ( tag       => 'value',
                                priority  => 'value',
                                annot     => 'value',
                                hugo      => 'value',
                                id        => 'value',
                                desc      => 'value',
                                hugo_non_noding => 'value',
                                id_non_coding   => 'value',
                                desc_non_coding => 'value',
                                hugo_3prim => 'value',
                                id_3prim   => 'value',
                                desc_3prim => 'value',
                                hugo_5prim => 'value',
                                id_5prim   => 'value',
                                desc_5prim => 'value',
                              );
  Exceptions  : none

=cut

sub getAnnotation{
  my $self = shift;
  my ($chr,$strand,$location,$sense,$tag) = @_;

  if(!defined $sense) {
    if($strand == 1) {
      $sense = 'after';
    } else {
      $sense = 'before';
    }
  }

  my $slice;

  eval{
    if (($sense eq 'before' && $strand == 1) || ($sense eq 'after' && $strand == -1)){
      $slice = $self->{slice_adaptor}->fetch_by_region( 'chromosome', $chr,$location-($self->{tag_length}-2), $location+1, $strand);
    }else{
      $slice = $self->{slice_adaptor}->fetch_by_region( 'chromosome', $chr,$location+1, $location+$self->{tag_length}, $strand);
    }
  }; 
  if ($@) {
    die ("annoteOneTag : $@");
    exit 1;
  }

  if (!$slice){
    die('annoteOneTag : connexion failed with tag at $location in chr$chromosome with strand = $srand_in');
    exit 1;
  }

  if(defined $tag && ($tag ne $slice->seq)) {
    print STDERR "original : $tag, ensembl : ". $slice->seq() . "\n";
    die('Sequence retrieved from Ensembl API is not consistant with the argument one.');
    exit 2;
  }

  return $self->annoteTagOnGenome($slice,$chr,$strand,$location);
}

=head1 PRIVATE METHODS

Please, do not try to use these methods outside this package.

=head2 annoteTagOnGenome

  Arg [1]     : Bio::EnsEMBL::Slice - Slice object use for annotation
  Arg [2]     : String - Chromosome
  Arg [3]     : String - Strand
  Arg [4]     : String - Position

  Exemple     : my $annotation = ReadAnnotation->annoteTagOnGenome();
  Description : Create an annotation hash for the given tag.
  ReturnType  : Annotation hash of a tag : 
                %annotation = ( tag       => 'value',
                                priority  => 'value',
                                annot     => 'value',
                                hugo      => 'value',
                                id        => 'value',
                                desc      => 'value',
                                hugo_non_noding => 'value',
                                id_non_coding   => 'value',
                                desc_non_coding => 'value',
                                hugo_3prim => 'value',
                                id_3prim   => 'value',
                                desc_3prim => 'value',
                                distance_of_3prim_gene => 'value',
                                hugo_5prim => 'value',
                                id_5prim   => 'value',
                                desc_5prim => 'value',
                                distance_of_5prim_gene => 'value',
                              );
  Exceptions  : none

=cut

sub annoteTagOnGenome {
  my $self = shift;
  my $slice = shift;
  my ($chr,$strand,$location) = @_;
  my %annotation;

  my $findNonCoding = 0;
  my (@genes,@transcripts,@exons);

  @genes = @{ $slice->get_all_Genes() };

  my $priority = CracTools::Config::getConfVar("MAX_PRIORITY");

  # Loop on genes
  foreach my $gene (@genes){
    if ($priority > 1){
      eval{
        @transcripts = @{ $gene->get_all_Transcripts() };
      };
      if ($@) {
        die ("annoteTagOnGenome : $@");
        exit 1;
      };    

      # Loop on transcripts
      foreach my $transcript (@transcripts){
        if ($priority > 1){
          eval{
            @exons = @{ $transcript->get_all_Exons() };
          }; 
          if ($@) {
            die ("annoteTagOnGenome : $@");
            exit 1;
          }   

          if (defined $gene->biotype()){
            #protein coding part
            if ($gene->biotype() eq "protein_coding" || $gene->biotype() eq "pseudogene" 
              || $gene->biotype =~ /IG_C/i || $gene->biotype =~ /IG_V/i 
              || $gene->biotype =~ /TR_V/i || $gene->biotype =~ /TR_C/i
              || $gene->biotype =~ /TR_J/i || $gene->biotype =~ /IG_D/i
              || $gene->biotype =~ /IG_J/i || $gene->biotype =~ /TR_D/i
              || $gene->biotype eq "polymorphic_pseudogene"){
              foreach my $exon (@exons){	
                if ($priority > 1){
                  my ($startP,$startCoding,$endP,$endCoding) = ($exon->start(),$exon->coding_region_start($transcript),$exon->end(),$exon->coding_region_end($transcript));
                  if ($gene->strand() == 1){
                    if ($startP <= 1 && ($endP-$self->{tag_length}) >= 0){
                      #				    if ($exon->end_phase() != -1){
                      # CDS_sense
                      if ($priority > 2 && defined $startCoding && $startCoding <= 1 && defined $endCoding && ($endCoding-$self->{tag_length}) >= 0){ 
                        $priority = 2;
                        $annotation{hugo} = $gene->external_name();
                        $annotation{desc} = $gene->biotype();
                        $annotation{id} = $transcript->stable_id();
                        # 5PRIM_UTR_sense
                      }elsif ($priority > 3 && defined $startCoding && $startCoding > 1){
                        $priority = 3;
                        $annotation{hugo} = $gene->external_name();
                        $annotation{desc} = $gene->biotype();
                        $annotation{id} = $transcript->stable_id();
# 3P                    RIM_UTR_sense	
                      }elsif ($priority > 1){
                        $priority = 1;
                        $annotation{hugo} = $gene->external_name();
                        $annotation{desc} = $gene->biotype();
                        $annotation{id} = $transcript->stable_id();
                      }
#INXON_sense    
                    }elsif ($priority > 4 
                      && ($endP > 1 && $startP < 1 && $endP <= $self->{tag_length} || $startP < $self->{tag_length} && $startP >= 1 && $endP > $self->{tag_length})){
                      $priority = 4;
                      $annotation{hugo} = $gene->external_name();
                      $annotation{desc} = $gene->biotype();
                      $annotation{id} = $transcript->stable_id();
#INTRON_sense
                    }elsif ($priority > 5){
                      $priority = 5;
                      $annotation{hugo} = $gene->external_name();
                      $annotation{desc} = $gene->biotype();
                      $annotation{id} = $transcript->stable_id();
                    }
                  }elsif ($priority > 6){
                    if ($startP <= 1 && ($endP-$self->{tag_length}) >= 0){
#				    if ($exon->end_phase() != -1){
# CDS_antisense
                      if ($priority > 7 && defined $startCoding && $startCoding <= 1 && defined $endCoding && ($endCoding-$self->{tag_length}) >= 0){ 
                        $priority = 7;
                        $annotation{hugo} = $gene->external_name();
                        $annotation{desc} = $gene->biotype();
                        $annotation{id} = $transcript->stable_id();
# 5PRIM_UTR_antisense	
                      }elsif ($priority > 8 && defined $startCoding && $startCoding > 1){
                        $priority = 8;		 
                        $annotation{hugo} = $gene->external_name();
                        $annotation{desc} = $gene->biotype();
                        $annotation{id} = $transcript->stable_id();
# 3PRIM_UTR_antisense
                      }elsif ($priority > 6){
                        $priority = 6;		 
                        $annotation{hugo} = $gene->external_name();
                        $annotation{desc} = $gene->biotype();
                        $annotation{id} = $transcript->stable_id();
                      }
# INXON_antisense
                    }elsif ($priority > 9
                      && ($endP > 1 && $startP < 1 && $endP <= $self->{tag_length} || $startP < $self->{tag_length} && $startP >= 1 && $endP > $self->{tag_length})){
                      $priority = 9;
                      $annotation{hugo} = $gene->external_name();
                      $annotation{desc} = $gene->biotype();
                      $annotation{id} = $transcript->stable_id();
# INTRON_antisense
                    }elsif ($priority > 10){
                      $priority = 10;
                      $annotation{hugo} = $gene->external_name();
                      $annotation{desc} = $gene->biotype();
                      $annotation{id} = $transcript->stable_id();
                    }
                  }
                }
              }
            }

#non-coding part
            if (!$findNonCoding && $gene->strand() == 1){
              if ($gene->biotype() eq "miRNA" || $gene->biotype() eq "miRNA_pseudogene" 
                || $gene->biotype() eq "snRNA" || $gene->biotype() eq "snRNA_pseudogene" 
                || $gene->biotype eq "snoRNA" || $gene->biotype() eq "snoRNA_pseudogene" 
                || $gene->biotype() eq "rRNA" || $gene->biotype() eq "rRNA_pseudogene"
                || $gene->biotype() eq "Mt_rRNA" || $gene->biotype() eq "Mt_rRNA_pseudogene"
                || $gene->biotype() eq "Mt_tRNA" || $gene->biotype() eq "Mt_tRNA_pseudogene"
                || $gene->biotype() eq "tRNA" || $gene->biotype() eq "tRNA_pseudogene"
                || $gene->biotype() eq "scRNA" || $gene->biotype() eq "scRNA_pseudogene"
                || $gene->biotype eq "ncRNA" || $gene->biotype eq "ncRNA_pseudogene"
                || $gene->biotype eq "3prime_overlapping_ncrna") {
                $annotation{desc_non_coding} = "small_ncRNA:".$gene->biotype();
                $findNonCoding = 1;
              }elsif ($gene->biotype() =~ /lincRNA/i){
                $annotation{desc_non_coding} = "lincRNA:".$gene->biotype();
                $findNonCoding = 1;
              }elsif ($gene->biotype() eq "antisense" || $gene->biotype() eq "sense_intronic"
                || $gene->biotype() eq "processed_transcript"){
                $annotation{desc_non_coding} = "other_lncRNA:".$gene->biotype();
                $findNonCoding = 1;
              }elsif ($gene->biotype() =~ /non_coding/ 
                || $gene->biotype() eq "misc_RNA" || $gene->biotype() eq "misc_RNA_pseudogene"
                || $gene->biotype() eq "ncrna_host" || $gene->biotype() eq "sense_overlapping"
                || $gene->biotype() eq "retained_intron"
                || $gene->biotype() eq "processed_pseudogene" || $gene->biotype() eq "unprocessed_pseudogene"
                || $gene->biotype() eq "transcribed_processed_pseudogene" 
                || $gene->biotype() eq "transcribed_unprocessed_pseudogene" 
                || $gene->biotype() eq "retrotransposed" || $gene->biotype() eq "unitary_pseudogene"
              ){
                $annotation{desc_non_coding} = "other_noncodingRNA:".$gene->biotype();
                $findNonCoding = 1;
              }elsif ($gene->biotype() ne "protein_coding" && $gene->biotype() ne "pseudogene" 
                && $gene->biotype !~ /IG_C/i && $gene->biotype !~ /IG_V/i 
                && $gene->biotype !~ /TR_V/i && $gene->biotype !~ /TR_C/i
                && $gene->biotype !~ /TR_J/i && $gene->biotype !~ /IG_D/i
                && $gene->biotype !~ /IG_J/i && $gene->biotype !~ /TR_D/i
                && $gene->biotype ne "polymorphic_pseudogene"){
                print STDERR "Warning: ".$gene->biotype." neither in process A nor process B!\n";  
              }
              if($findNonCoding) {
                $annotation{hugo_non_coding} = $gene->external_name();
                $annotation{id_non_coding} = $transcript->stable_id();
              }
            }
          }else{
            print STDERR "Warning: no biotype!\n";  
          }
        }
      }
    }
  }

  # If we haven'y fin anything intersting, we look a little bit further
  if ($priority == CracTools::Config::getConfVar("MAX_PRIORITY")){
    my ($sliceEST, $sliceProxy5PRIM, $sliceProxy3PRIM);
    eval{
      if ($strand == -1){
        $sliceProxy5PRIM = $self->{slice_adaptor}->fetch_by_region( 'chromosome', $chr,$location-($self->{tag_length}-1), $location + $self->{DISTANCE_MAX}, $strand);
        $sliceProxy3PRIM = $self->{slice_adaptor}->fetch_by_region( 'chromosome', $chr,$location - $self->{DISTANCE_MAX}, $location-($self->{tag_length}-1),  $strand);
        #$sliceEST = $self->{sliceEST_adaptor}->fetch_by_region( 'chromosome', $chr,$location-($self->{tag_length}-1), $location, $strand);
      }else{
        $sliceProxy5PRIM = $self->{slice_adaptor}->fetch_by_region( 'chromosome', $chr,$location + 1 - $self->{DISTANCE_MAX}, $location+$self->{tag_length}, $strand);
        $sliceProxy3PRIM = $self->{slice_adaptor}->fetch_by_region( 'chromosome', $chr, $location+$self->{tag_length}, $location + 1 + $self->{DISTANCE_MAX}, $strand);
        #$sliceEST = $self->{sliceEST_adaptor}->fetch_by_region( 'chromosome', $chr,$location+1, $location+$self->{tag_length}, $strand);
      }
    };
    if ($@) {
      die ("printAnnotation : $@");
      exit 1;
    }    

    my ($findEST,$findProxy) = (0,0);
    my $posGene = $self->{DISTANCE_MAX};
    my $posGeneTmp;

    #if ($#{ $sliceEST->get_all_Genes() } != -1){
    #  @genes = @{ $sliceEST->get_all_Genes() };
    #  $annotation{id} = $genes[0]->stable_id();
    #  $findEST = 1;
    #}

    if($#{ $sliceProxy5PRIM->get_all_Genes() } != -1){
      @genes = @{ $sliceProxy5PRIM->get_all_Genes() };
      foreach my $gene (@genes){
        if ($gene->strand() == 1){
          $posGeneTmp = $self->{DISTANCE_MAX} - $gene->end();
          if ($posGeneTmp < $posGene && ($self->{DISTANCE_MAX} - $gene->end()) > 0){
            $annotation{hugo_5prim} = $gene->external_name();
            $annotation{desc_5prim} = $gene->biotype();
            $annotation{id_5prim} = $gene->stable_id();
            $posGene = $posGeneTmp;
            $annotation{distance_of_5prim_gene} = $posGene;
            $findProxy = 1;
          }
        }
      }
    }

    if ($findProxy && $posGene <= $self->{INTERGENIC_THRESHOLD}){
      $priority = 11;
    } else{
      if ($findEST){
        $priority = 12;
      } else{
        $priority = 13;
      }
    }

    $findProxy = 0;
    if($#{ $sliceProxy3PRIM->get_all_Genes() } != -1){
      @genes = @{ $sliceProxy3PRIM->get_all_Genes() };
      foreach my $gene (@genes){
        if ($gene->strand() == 1){
          $posGeneTmp = $gene->start();
          if ($posGeneTmp < $posGene && $gene->start() > 0){
            $annotation{hugo_3prim} = $gene->external_name();
            $annotation{desc_3prim} = $gene->biotype();
            $annotation{id_3prim} = $gene->stable_id();
            $posGene = $posGeneTmp;
            $annotation{distance_of_3prim_gene} = $posGene;
            $findProxy = 1;
          }	
        }
      }
    }
  }

  $annotation{priority} = $priority;
  $annotation{tag} = $slice->seq();
  $annotation{annot} = $priority_annotation{$priority};

  return \%annotation;
}

1;
