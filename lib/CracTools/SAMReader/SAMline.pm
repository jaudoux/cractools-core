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

=encoding utf8

=head1 NAME

CracTools::SAMReader::SAMline - The object for manipulation a SAM line.

=head1 SYNOPSIS

  use CracTools::SAMReader::SAMline;

  $sam_line = CracTools::SAMReader::SAMline->new($line);

=head1 DESCRIPTION

An object for easy acces to SAM line fields. See SAM Specifications for more informations :
http://samtools.sourceforge.net/SAM1.pdf

=cut

package CracTools::SAMReader::SAMline;


use strict;
use warnings;
use Carp;
use Data::Dumper;
use Switch;

use CracTools::Utils;

=head1 Variables

=head2 %flags

SAM flags :

=over 2

=item * MULTIPLE_SEGMENTS => 1

=item * PROPERLY_ALIGNED => 2

=item * UNMAPPED => 4,

=item * NEXT_UNMAPPED => 8,

=item * REVERSE_COMPLEMENTED => 16,

=item * NEXT_REVERSE_COMPLEMENTED => 32,

=item * FIRST_SEGMENT => 64,

=item * LAST_SEGMENT => 128,

=item * SECONDARY_ALIGNMENT => 256,

=item * QUALITY_CONTROLS_FAILED => 512,

=item * PCR_DUPLICATED => 1024,

=item * CHIMERIC_ALIGNMENT => 2048,

=back

=cut

our %flags = ( MULTIPLE_SEGMENTS => 1,
            PROPERLY_ALIGNED => 2,
            UNMAPPED => 4,
            NEXT_UNMAPPED => 8,
            REVERSE_COMPLEMENTED => 16,
            NEXT_REVERSE_COMPLEMENTED => 32,
            FIRST_SEGMENT => 64,
            LAST_SEGMENT => 128,
            SECONDARY_ALIGNMENT => 256,
            QUALITY_CONTROLS_FAILED => 512,
            PCR_DUPLICATED => 1024,
            CHIMERIC_ALIGNMENT => 2048,
          );

=head1 STATIC PARSING METHODS

These methods can be used without creating an CracTools::SAMReader::SAMline object.
They are designed to provided efficient performance when parsing huge SAM files,
because creating object in Perl can be long and useless for some purposes.

=cut

=head2 hasEvent

  Arg [1] : String - SAM line
  Arg [2] : eventType

=cut

sub hasEvent {
  my ($line,$event_type) = @_;
  croak("Missing argument(s)") unless defined $line && defined $event_type;
  return $line =~ /XE:Z:\d+:\d+:$event_type/i;
}

=head1 Methods

=head2 new

  Arg [1] : String - SAM line in TAB-separated format.

  Example     : $sam_line = CracTools::SAMline->new$($line);
  Description : Create a new CracTools::SAMline obect.
  ReturnType  : CracTools::SAMline
  Exceptions  : none

=cut

sub new {
  my $class = shift;
  my ($line) = @_;
  chomp $line;  
  my ($qname,$flag,$rname,$pos,$mapq,$cigar,$rnext,$pnext,$tlen,$seq,$qual,@others) = split("\t",$line);
  my %extended_fields;
  foreach(@others) {
    my $key = substr($_,0,2);
    my $value = substr($_,5);
    #my($key,$type,$value) = split(':',$_,3); 
    # Hack in case there is multiple field with the same tag
    $extended_fields{$key} = defined $extended_fields{$key}? $extended_fields{$key}.";".$value : $value;
  }

  my $self = bless{ 
    qname => $qname,
    flag => $flag,
    rname => $rname,
    pos => $pos,
    mapq => $mapq,
    cigar => $cigar,
    rnext => $rnext,
    pnext => $pnext,
    tlen => $tlen,
    seq => $seq,
    qual => $qual,
    extended_fields => \%extended_fields,
    line => $line,
  }, $class;

  return $self;
}

=head2 isFlagged

  Arg [1] : Integer - The flag to test (1,2,4,8, ... ,1024)

  Example     : if($SAMline->isFlagged($fags{unmapped}) {
                  DO_SOMETHING... 
                };
  Description : Test if the line has the flag in parameter setted.
  ReturnType  : Boolean
  Exceptions  : none

=cut

sub isFlagged {
  my $self = shift;
  my $flag = shift;
  return $self->flag & $flag;
}

=head2 getStrand

  Example     : $strand = $SAMline->getStrand(); 
  Description : Return the strand of the SAMline :
                - "1" if forward strand
                - "-1" if reverse strand
  ReturnType  : 1 or -1
  Exceptions  : none

=cut

sub getStrand {
  my $self = shift;
  if($self->isFlagged($flags{REVERSE_COMPLEMENTED})) {
    return -1;
  } else {
    return 1;
  }
}

=head2 getOriginalSeq

  Descrition   : Return the original sequence as it was in the FASTQ file.
                 In fact we reverse complemente the sequence if flag 16 is raised.

=cut

sub getOriginalSeq {
  my $self = shift;
  if($self->isFlagged($flags{REVERSE_COMPLEMENTED})) {
    return CracTools::Utils::reverseComplement($self->seq);
  } else {
    return $self->seq;
  }
}

=head2 getLocAsCracFormat

  Example     : $loc = $SAMline->getLocAsCracFormat(); 
  Description : Return the location of the sequence using CRAC format : "chr|strand,position".
                For example : X|-1,2154520
  ReturnType  : String
  Exceptions  : none

=cut

sub getLocAsCracFormat {
  my $self = shift;
  return $self->rname."|".$self->getStrand.",".$self->pos;
}

=head2 getPatch

  Description : If the SAMline has been modified, this method will generate
                a patch in UnifiedDiff format that represent the changes.
  ReturnType  : String (patch) if line has changed, False (0) either.
  Exceptions  : none

=cut

sub getPatch {
  my $self = shift;
  my $line_number = shift;
  croak("Cannot generate patch without the line number in argument") unless defined $line_number;
  if($self->updatedLine ne $self->line) {
    my $line1 = $self->line;
    my $line2 = $self->updatedLine;
    chomp($line1);
    chomp($line2);
    return "@@ -$line_number,1 +$line_number,1 @@\n-$line1\n+$line2";
  } else {
    return 0;
  }
}

=head1 GETTERS AND SETTERS

=cut

=head2 line

  Description : Getterr for the whole SAMline as a string.
  ReturnType  : String
  Exceptions  : none

=cut

sub line {
  my $self = shift;
  return $self->{line};
}

=head2 updatedLine

  Description : Getter/Setter for the updated line.
                If there is not updated line, this method return
                the original SAM line.
  RetrunType  : String

=cut

sub updatedLine {
  my $self = shift;
  my $updated_line = shift;
  if(defined $updated_line) {
    $self->{updated_line} = $updated_line;
    return 1;
  } elsif (!defined $self->{updated_line}) {
    return $self->line;
  } else {
    return $self->{updated_line};
  }
}

=head2 qname

  Description : Getter/Setter for attribute qname
  ReturnType  : String
  Exceptions  : none

=cut

sub qname {
  my $self = shift;
  my $new_qname = shift;
  if(defined $new_qname) {
    $self->{qname} = $new_qname;
  }
  return $self->{qname};
}

=head2 flag

  Description : Getter/Setter for attribute flag
  ReturnType  : String
  Exceptions  : none

=cut

sub flag {
  my $self = shift;
  my $new_flag = shift;
  if(defined $new_flag) {
    $self->{flag} = $new_flag;
  }
  return $self->{flag};
}

=head2 rname

  Description : Getter/Setter for attribute rname (chromosome for eucaryotes)
  ReturnType  : String
  Exceptions  : none

=cut

sub rname {
  my $self = shift;
  my $new_rname = shift;
  if(defined $new_rname) {
    $self->{rname} = $new_rname;
  }
  return $self->{rname};
}

=head2 chr

  Description : Getter/Setter for attribute rname (Alias)
  ReturnType  : String
  Exceptions  : none

=cut

sub chr {
  my $self = shift;
  $self->rname(@_);
}

=head2 pos

  Description : Getter/Setter for attribute pos (position of the sequence)
  ReturnType  : String
  Exceptions  : none

=cut

sub pos {
  my $self = shift;
  my $new_pos = shift;
  if(defined $new_pos) {
    $self->{pos} = $new_pos;
  }
  return $self->{pos};
}

=head2 mapq

  Description : Getter/Setter for attribute mapq (mapping quality)
  ReturnType  : String
  Exceptions  : none

=cut

sub mapq {
  my $self = shift;
  my $new_mapq = shift;
  if(defined $new_mapq) {
    $self->{mapq} = $new_mapq;
  }
  return $self->{mapq};
}

=head2 cigar

  Description : Getter/Setter for attribute cigar (see SAM doc)
  ReturnType  : String
  Exceptions  : none

=cut

sub cigar {
  my $self = shift;
  my $new_cigar = shift;
  if(defined $new_cigar) {
    $self->{cigar} = $new_cigar;
  }
  return $self->{cigar};
}

=head2 rnext

  Description : Getter/Setter for attribute rnext (see SAM doc)
  ReturnType  : String
  Exceptions  : none

=cut

sub rnext {
  my $self = shift;
  my $new_rnext = shift;
  if(defined $new_rnext) {
    $self->{rnext} = $new_rnext;
  }
  return $self->{rnext};
}

=head2 pnext

  Description : Getter/Setter for attribute pnext (see SAM doc)
  ReturnType  : Integer
  Exceptions  : none

=cut

sub pnext {
  my $self = shift;
  my $new_pnext = shift;
  if(defined $new_pnext) {
    $self->{pnext} = $new_pnext;
  }
  return $self->{pnext};
}

=head2 tlen

  Description : Getter/Setter for attribute tlen (sequence length)
  ReturnType  : Integer
  Exceptions  : none

=cut

sub tlen {
  my $self = shift;
  my $new_tlen = shift;
  if(defined $new_tlen) {
    $self->{tlen} = $new_tlen;
  }
  return $self->{tlen};
}

=head2 seq

  Description : Getter/Setter for attribute seq (the sequence).
                Please use getOriginalSeq if you want to retrieve the oriented
                sequence, that what you need in most cases.
  ReturnType  : String
  Exceptions  : none

=cut

sub seq {
  my $self = shift;
  my $new_seq = shift;
  if(defined $new_seq) {
    $self->{seq} = $new_seq;
  }
  return $self->{seq};
}

=head2 qual

  Description : Getter/Setter for attribute qual (sequence quality)
  ReturnType  : String
  Exceptions  : none

=cut

sub qual {
  my $self = shift;
  my $new_qual = shift;
  if(defined $new_qual) {
    $self->{qual} = $new_qual;
  }
  return $self->{qual};
}

=head2 getOptionalField

  Example     : 
  Description : 
  ReturnType  : 

=cut

sub getOptionalField {
  my $self = shift;
  my $field = shift;
  croak("Missing \$field argument to call getOptionalField") unless defined $field;
  return $self->{extended_fields}{$field};
}


=head2 getChimericAlignments

  Description : Parser of SA fields of SAM file in order to find chimeric reads
  ReturnType  : Array reference
                Elements are hash [ chr    => String, 
                                    pos    => int, 
                                    strand => 1/-1, 
                                    cigar  => String,
                                    mapq   => int,
                                    edist  => int
                                  ]

=cut

sub getChimericAlignments {
    my $self = shift;
    # check the existence of the SA field in the SAM line
    if (defined $self->{extended_fields}{SA}){
	my @array_hash;
	my (@SA_alignments) = split(/;/,$self->{extended_fields}{SA});
	for (my $i=0 ;  $i < scalar @SA_alignments ; $i++){
	    my ($chr,$pos,$strand,$cigar,$mapq,$edist) = split(/,/,$SA_alignments[$i]);
	    # strand switch from "+,-" to "1,-1"
	    if ($strand eq '+'){
		$strand = 1;
	    }else{
		$strand = -1;
	    }
	    my $hash = { chr => $chr, 
			 pos => $pos,
			 strand => $strand,
			 cigar => $cigar,
			 mapq => $mapq,
			 edist => $edist};
	    push(@array_hash,$hash);
	}
	return \@array_hash;
    }
    return undef;
}



=head2 pSupport

  Description : Return the support profile of the read if the SAM file has been generated with
                CRAC option --detailed
  ReturnType  : String

=cut

sub pSupport {
  my $self = shift;
  $self->loadSamDetailed;
  return $self->{sam_detailed}{p_support};
}

=head2 pLoc

  Description : Return the location profile of the read if the SAM file has been generated with
                CRAC option --detailed
  ReturnType  : String

=cut

sub pLoc {
  my $self = shift;
  $self->loadSamDetailed;
  return $self->{sam_detailed}{p_loc};
}

=head2 pairedChimera

  Description : return the chimeric coordinates of the paired chimera associated to this read if there is one

  ReturnType  : array(chr1,pos1,strand1,chr2,pos2,strand2) or undef

=cut

sub pairedChimera {
  my $self = shift;
  $self->loadPaired();
  if(defined $self->{extended_fields}{XP}{chimera}) {
    my ($crac_loc1,$crac_loc2) = split(":",$self->{extended_fields}{XP}{chimera});
    return (expandCracLoc($crac_loc1),expandCracLoc($crac_loc2));
  } else {
    return undef;
  }
}

=head2 isPairedClassified

  Arg [1] : String - The class to test :
            - "unique"
            - "duplicated"
            - "multiple"

  Description : Test paired-end read clasification

  ReturnType  : Boolean

=cut

sub isPairedClassified {
  my $self = shift;
  my $class = shift;
  $self->loadPaired();

  if(defined $self->{extended_fields}{XP}{loc} && ref($self->{extended_fields}{XP}{loc}) ne 'HASH') {
    my ($uniq,$dupli,$multi) = split(":",$self->{extended_fields}{XP}{loc});
    $self->{extended_fields}{XP}{loc} = {unique => $uniq, duplicated => $dupli, multiple => $multi};

  }
  return $self->{extended_fields}{XP}{loc}{$class};
}


=head2 genericInfo
  
  [1] : Key of the generic info
  [2] : (Optional) Value of the generic info

  Description : Getter/Setter enable to store additional (generic) information 
                about the SAMline as a Key/Value. 
  Example : # Set a generic info
            $read->genericInfo("foo","bar")

            # Get a generic info
            print $read->genericInfo("foo"); # this will print "bar"
  ReturnType : ?
  Exceptions : none

=cut

sub genericInfo {
  my ($self,$key,$value) = @_;
  if(!defined $key) {
    die "You need to provide a key in order to set or get a genericInfo\n";
  } elsif(defined $value) {
    $self->{genericInfo}{$key} = $value;
  } else {
    return $self->{genericInfo}{$key};
  }
}

=head2 isClassified

  Arg [1] : String - The class to test :
            - "unique"
            - "duplicated"
            - "multiple"
            - "normal"
            - "almostNormal"

  Example     : if($sam_line->isClassified('normal')) {
                  DO_SOMETHING;
                }
  Description : Test if the line is classified according to the parameter value.
  ReturnType  : Boolean
  Exceptions  : none

=cut

sub isClassified {
  my $self = shift;
  my $class = shift;

  switch($class) {
    case "unique"       { return $self->{extended_fields}{XU} }
    case "duplicated"   { return $self->{extended_fields}{XD} }
    case "multiple"     { return $self->{extended_fields}{XM} }
    case "normal"       { return $self->{extended_fields}{XN} == 1 }
    case "almostNormal" { return $self->{extended_fields}{XU} == 2 }
    else                { return undef }
  }
}

=head2 events

  Arg [1] : String - The event type to return :
            - Junction
            - Ins
            - Del
            - SNP
            - Error
            - Chimera
            - Undetermined
            - BioUndetermined
            - ... (see CRAC SAM format specifications for more informations).
  Example     : my @junctions = @{$line->events('Junction')};
                foreach my $junction (@junctions) {
                  print "Foud Junction : [type : $junction->{type}, loc : $junction->{loc}, gap : $junction->{gap}]\n";
                } 
  Description : Return all events of the type specified in parameter
  ReturnType  : Array reference
  Exceptions  : none

=cut

sub events {
  my $self = shift;
  my $event_type = shift;
  $self->loadEvents($event_type);
  if(defined $self->{events}{$event_type}) {
    return $self->{events}{$event_type};
  } else {
    return [];
  }
}

=head1 PRIVATE METHODS

=head2 loadEvents

  Example     : $sam_line->loadEvents();
  Description : Loading of events attributes
  ReturnType  : none
  Exceptions  : none

=cut

sub loadEvents {
  my $self = shift;
  my $event_type_to_load = shift;
  ## TODO avoid double loading when doing lazy loading
  if(defined $event_type_to_load && defined $self->{$event_type_to_load}{loaded}) {
    return 0;
  }
  if(!defined $self->{events} && defined $self->{extended_fields}{XE}) { 
    # Init events
    my @events = split(";",$self->{extended_fields}{XE});
    foreach my $event (@events) {
      my ($event_id,$event_break_id,$event_type,$event_infos) = $event =~ /([^:]+):([^:]+):([^:]+):(.*)/g;
      next if(defined $event_type_to_load && $event_type ne $event_type_to_load);
      if(defined $event_id) {
        my %event_hash;
        if($event_type eq 'Junction') {
          my ($type,$pos_read,$loc,$gap) = split(':',$event_infos);
          my ($chr,$pos,$strand) = expandCracLoc($loc);
          %event_hash = ( type => $type,
                           pos => $pos_read,
                           loc => {chr => $chr, pos => $pos, strand => $strand},
                           gap => $gap,
                         );
        } elsif($event_type eq 'Ins' || $event_type eq 'Del') {
          my ($score,$pos_read,$loc,$nb) = split(':',$event_infos); 
          my ($chr,$pos,$strand) = expandCracLoc($loc);
          %event_hash = ( score => $score,
                      pos => $pos_read,
                      loc => {chr => $chr, pos => $pos, strand => $strand},
                      nb => $nb,
                    );
        } elsif($event_type eq 'SNP') {
          my ($score,$pos_read,$loc,$expected,$actual) = split(':',$event_infos); 
          my ($chr,$pos,$strand) = expandCracLoc($loc);
          %event_hash = ( score => $score,
                      pos => $pos_read,
                      loc => {chr => $chr, pos => $pos, strand => $strand},
                      expected => $expected,
                      actual => $actual,
                    );
        } elsif($event_type eq 'Error') {
          my ($type,$pos,$score,$other1,$other2) = split(':',$event_infos); 
          %event_hash = ( score => $score,
                      pos => $pos,
                      type => $type,
                      other1 => $other1,
                      other2 => $other2,
                    );
        } elsif($event_type eq 'chimera') {
          my ($pos_read,$loc1,$loc2) = split(':',$event_infos); 
          my ($chr1,$pos1,$strand1) = expandCracLoc($loc1);
          my ($chr2,$pos2,$strand2) = expandCracLoc($loc2);
          %event_hash = ( pos => $pos_read,
                      loc1 => {chr => $chr1, pos => $pos1, strand => $strand1},
                      loc2 => {chr => $chr2, pos => $pos2, strand => $strand2},
                    );
        } elsif($event_type eq 'Undetermined') {
          %event_hash = ( message => $event_infos,
                    );
        } elsif($event_type eq 'BioUndetermined') {
          my ($pos,$message) = $event_infos =~ /([^:]+):(.*)/; 
          %event_hash = ( pos => $pos,
                      message => $message,
                    );
        }
        if (keys %event_hash > 1) {
          $event_hash{event_id} = $event_id;
          $event_hash{break_id} = $event_break_id;
          $event_hash{event_type} = $event_type;
          $self->addEvent(\%event_hash);
        }
      }
    }
    # If we have only load a specific event type
    if(defined $event_type_to_load) {
      $self->{$event_type_to_load}{loaded} = 1;
    # Else we have load every events.
    } else {
      $self->{events}{loaded} = 1;
    }
  }
}

=head2 addEvent

  Arg [1] : String - The event type
  Arg [2] : Hash reference - The event object
  Example     : $line->addEvent($event_type,\%event); 
  Description : Return all events of the type specified in parameter
  ReturnType  : none
  Exceptions  : none

=cut

sub addEvent {
  my $self = shift;
  my $event = shift;
  my $event_type = $event->{event_type};
  if(defined $self->{events}{$event_type}) {
    push(@{$self->{events}{$event_type}},$event);
  } else {
    $self->{events}{$event_type} = [$event];
  }
}

=head2 removeEvent 

  Arg [1] : Hash reference - The event object

  Description : Remove the event from the event hash and from the line.

=cut

sub removeEvent {
  my $self = shift;
  my $delete_event = shift;
  my $type = $delete_event->{event_type};
  if(defined $type && defined $self->{events}{$type}) {
    my $i = 0;
    foreach my $event (@{$self->{events}{$type}}) {
      if($event eq $delete_event) {
        splice @{$self->{events}{$type}}, $i, 1;
        return 1;
      }
      $i++;
    }
  }
  return 0;
}

=head2 updateEvent
  

=cut

sub updateEvent {
  my $self = shift;
  my $event = shift;
  my $new_event_type = shift;
  my %new_event = @_;

  # Update new event with old break id and event id
  $new_event{event_type} = $new_event_type;
  $new_event{event_id} = $event->{event_id};
  $new_event{break_id} = $event->{break_id};

  if($self->removeEvent($event)) {
    # Catch warnings on string concatenation that correspond to missing
    # field in the hash for the event to update
    local $SIG{'__WARN__'} = sub {croak("Invalid event hash for event type '$new_event_type'");};
    my $base_XE = 'XE:Z:'.$new_event{event_id}.':'.$new_event{break_id};
    my $new_XE = $base_XE . ':';
    if($new_event_type eq 'Junction') {
        my $loc = compressCracLoc($new_event{loc}{chr},$new_event{loc}{pos},$new_event{loc}{strand});
        $new_XE .= $new_event_type.':'.
                   $new_event{type}.':'.
                   $new_event{pos}.':'.
                   $loc.':'.
                   $new_event{gap};
    } elsif($new_event_type eq 'Ins' || $new_event_type eq 'Del') {
        my $loc = compressCracLoc($new_event{loc}{chr},$new_event{loc}{pos},$new_event{loc}{strand});
        $new_XE .= $new_event_type.':'.
                   $new_event{score}.':'.
                   $new_event{pos}.':'.
                   $loc.':'.
                   $new_event{nb};
    } elsif($new_event_type eq 'SNP') {
        my $loc = compressCracLoc($new_event{loc}{chr},$new_event{loc}{pos},$new_event{loc}{strand});
        $new_XE .= $new_event_type.':'.
                   $new_event{score}.':'.
                   $new_event{pos}.':'.
                   $loc.':'.
                   $new_event{expected}.':'.
                   $new_event{actual};
    } elsif($new_event_type eq 'Error') {
        my $loc = compressCracLoc($new_event{loc}{chr},$new_event{loc}{pos},$new_event{loc}{strand});
        $new_XE .= $new_event_type.':'.
                   $new_event{type}.':'.
                   $new_event{pos}.':'.
                   $new_event{score}.':'.
                   $new_event{other1}.':'.
                   $new_event{other2};
    } elsif($new_event_type eq 'chimera') {
        my $loc1 = compressCracLoc($new_event{loc1}{chr},$new_event{loc1}{pos},$new_event{loc1}{strand});
        my $loc2 = compressCracLoc($new_event{loc2}{chr},$new_event{loc2}{pos},$new_event{loc2}{strand});
        $new_XE .= $new_event_type.':'.
                   $new_event{pos}.':'.
                   $loc1.':'.
                   $loc2;
    } elsif($new_event_type eq 'Undetermined') {
        $new_XE .= $new_event_type.':'.
                   $new_event{message};
    } elsif($new_event_type eq 'BioUndetermined') {
        $new_XE .= $new_event_type.':'.
                   $new_event{pos}.':'.
                   $new_event{message};
    } else {
      croak("Unknown type of event : $new_event_type");
    }
    $self->addEvent(\%new_event);
    my $new_line = $self->updatedLine;
    $new_line =~ s/($base_XE:[^\t]*)/$new_XE/;
    $self->updatedLine($new_line);
  } else {
    croak('Event not find');
  }
}

=head2 loadSamDetailed

  Example     : $sam_line->loadSamDetailed();
  Description : Loading of sam detaileds attributes
  ReturnType  : none
  Exceptions  : none

=cut

sub loadSamDetailed {
  my $self = shift;
  if(!defined $self->{sam_detailed}) {
    #my ($detailed) = $self->line =~ /XR:Z:([^\s]+)/g;
    my $detailed = $self->{extended_fields}{XR};
    my @detailed_fields = split(";",$detailed); 
    foreach (@detailed_fields) {
      my ($k,$v) = split('=',$_);
      if($k eq 'p_loc') {
        $self->{sam_detailed}{p_loc} = $v;
      } elsif($k eq 'p_support') {
        $self->{sam_detailed}{p_support} = $v;
      } else {
        carp("Unknown sam detailed field : $k");
      }
    }
    $self->{sam_detailed}{loaded} = 1;
  }
}

=head2 loadPaired

  Example     : $sam_line->loadPaired();
  Description : Loading of sam detaileds attributes
  ReturnType  : none
  Exceptions  : none

=cut

sub loadPaired {
  my $self = shift;
  # If XP fields exist and we haven't load it already
  if(defined $self->{extended_fields}{XP} && ref($self->{extended_fields}{XP}) ne 'HASH') {
    my @paired_fields = split(";",$self->{extended_fields}{XP});
    $self->{extended_fields}{XP} = {}; # Chamge type of XP exetended field from scalar to hash ref
    foreach (@paired_fields) {
      my($key,$value) = split(":",$_,2);
      $self->{extended_fields}{XP}{$key} = $value;
    }
  }
}

=head2 expandCracLoc

  Arg [1] : String - Localisation in crac format : Chromosome|strand,position
            Ex : X|-1,2332377

  Description : Extract Chromosme, position and strand as separated variable from
                the localisation in CRAC format.
  ReturnType  : Array($chromosome,$position,$strand)

=cut

sub expandCracLoc {
  my $loc = shift;
  my($chr,$strand,$pos) = $loc =~ /(\S+)\|(\S+)?,(\S+)?/; 
  return ($chr,$pos,$strand);
}

=head2 compressCracLoc

  Arg [1] : String - Chromosome
  Arg [2] : Integer - Postition
  Arg [3] : Integer (1,-1) - Strand

  Description : Reverse function of "expandCracLoc"
  ReturnType  : String (localisation in CRAC format)

=cut

sub compressCracLoc {
  my ($chr,$pos,$strand) = @_;
  confess("Missing argument") unless defined $chr && defined $pos && defined $strand;
  return $chr."|".$strand.",".$pos;
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
