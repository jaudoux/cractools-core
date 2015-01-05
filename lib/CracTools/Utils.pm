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
use Fcntl qw( SEEK_SET );

=head1 SYNOPSIS

  # Reverse complementing a sequence
  my $seq = reverseComplemente("ATGC");

  # Reading a FASTQ file
  my $it = seqFileIterator('file.fastq','fastq');
  while(my $entry = $it->()) {
    print "Sequence name   : $entry->{name}
           Sequence        : $entry->{seq}
           Sequence quality: $entry->{qual}","\n";
  }

  # Reading paired-end files easier
  my $it = pairedEndSeqFileIterator($file);
  while (my $entry = $it->()) {
    print "Read_1 : $entry->{read1}->{seq}
           Read_2 : $entry->{read2}->{seq}";
  }

  # Parsing a GFF file
  my $it = gffFileIterator($file);
  while (my $annot = $it->()) {
    print "chr    : $annot->{chr}
           start  : $annot->{start}
           end    : $annot->{end}";
  }


=head1 DESCRIPTION

Bio::Lite is a set of subroutines that aims to answer similar questions as
Bio-perl distribution in a FAST and SIMPLE way.

Bio::Lite does not make use of complexe data struture, or
objects, that would lead to a slow execution.

All methods can be imported with a single "use Bio::Lite".

Bio::Lite is a lightweight-single-module with NO DEPENDENCIES.

=head1 UTILS

=head2 reverseComplement

Reverse complemente the (nucleotid) sequence in arguement.

Example:
  
  my $seq_revcomp = reverseComplement($seq);

reverseComplement is more than B<100x faster than Bio-Perl> revcom_as_string()

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

Convert strand from '+/-' standard to '1/-1' standard and the opposite.

Example:

  say "Forward a: ",convertStrand('+');
  say "Forward b: ",convertStrand(1);
  say "Reverse a: ",convertStrand('-');
  say "Reverss b: ",convertStrand(-1);

will print

  Forward a: 1
  Forward b: +
  Reverse a: -1
  Reverse b: -

=cut

my %conversion_hash = ( '+' => 1, '-' => '-1', '1' => '+', '-1' => '-');
sub convertStrand($) {
  my $strand = shift;
  return $conversion_hash{$strand};
}

=head1 PARSING

This are some tools that aim to read (bio) files like

=over

=item Sequence files : FASTA, FASTQ

=item Annotation files : GFF3, GTF2, BED6, BED12, ...

=item Alignement files : SAM, BAM

=back

=head2 seqFileIterator

Open Fasta, or Fastq files (can be gziped).
seqFileIterator has an automatic file extension detection but you can force it
using a second parameter with the format : 'fasta' or 'fastq'.

Example:

  my $it = seqFileIterator('file.fastq','fastq');
  while(my $entry = $it->()) {
    print "Sequence name   : $entry->{name}
           Sequence        : $entry->{seq}
           Sequence quality: $entry->{qual}","\n";
  }

Return: HashRef
  
  { name => 'sequence_identifier',
    seq  => 'sequence_value',
    qual => 'sequence_quality', # only defined for FASTQ files
  }

seqFileIterator is more than B<50x faster than Bio-Perl> Bio::SeqIO for FASTQ files
seqFileIterator is 4x faster than Bio-Perl Bio::SeqIO for FASTA files

=cut

sub seqFileIterator {
  my ($file,$format) = @_;

  croak "Missing file in argument of seqFileIterator" if !defined $file;

  # Get file handle for $file
  my $fh = getReadingFileHandle($file);

  # Automatic file extension detection
  if(!defined $format) {
    if($file =~ /\.(fasta|fa)(\.|$)/) {
      $format = 'fasta';
    } elsif($file =~ /\.(fastq|fq)(\.|$)/) {
      $format = 'fastq';
    } else {
      croak "Undefined file extension";
    }
  } else {
    $format = lc $format;
  }

  # FASTA ITERATOR
  if ($format eq 'fasta') {
    # Read prev line for FASTA because we dont know the number
    # of line used for the sequence
    my $prev_line = <$fh>;
    chomp $prev_line;
    return sub {
      my ($name,$seq,$qual); 
      if(defined $prev_line) {
        ($name) = $prev_line =~ />(.*)$/;
        $prev_line = <$fh>;
        # Until we find a new sequence identifier ">", we
        # concatenate the lines corresponding to the sequence
        while(defined $prev_line && $prev_line !~ /^>/) {
          chomp $prev_line;
          $seq .= $prev_line;
          $prev_line = <$fh>;
        }
        return {name => $name, seq => $seq, qual => $qual};
      } else {
        return undef;
      }
    };
  # FASTQ ITERATOR
  } elsif ($format eq 'fastq') {
    return sub {
      my ($name,$seq,$qual); 
      my $line = <$fh>;
      ($name) = $line =~ /@(.*)$/ if defined $line;
      if(defined $name) {
        $seq = <$fh>;
        chomp $seq;
        <$fh>; # skip second seq name (useless line)
        $qual = <$fh>;
        chomp $qual;
        return {name => $name, seq => $seq, qual => $qual};
      } else {
        return undef;
      }
    };
  } else {
    croak "Undefined file format";
  }
}

=head2 pairedEndSeqFileIterator

Open Paired-End Sequence files using seqFileIterator()

Paird-End files are generated by Next Generation Sequencing technologies (like Illumina) where two
reads are sequenced from the same DNA fragment and saved in separated files.

Example:

  my $it = pairedEndSeqFileIterator($file);
  while (my $entry = $it->()) {
    print "Read_1 : $entry->{read1}->{seq}
           Read_2 : $entry->{read2}->{seq}";
  }

Return: HashRef

  { read1 => 'see seqFileIterator() return',
    read2 => 'see seqFileIterator() return'
  }

pairedEndSeqFileIterator has no equivalent in Bio-Perl

=cut

sub pairedEndSeqFileIterator {
  my($file1,$file2,$format) = @_;

  my $it1 = seqFileIterator($file1,$format);
  my $it2 = seqFileIterator($file2,$format);

  return sub {
    my $entry1 = $it1->();
    my $entry2 = $it2->();
    if(defined $entry1 && defined $entry2) {
      return { read1 => $entry1, read2 => $entry2 };
    } else {
      return undef;
    }
  };
}

=head2 bedFileIterator 

manage BED files format

Example:

  my $it = bedFileIterator($file);
  while (my $annot = $it->()) {
    print "chr    : $annot->{chr}
           start  : $annot->{start}
           end    : $annot->{end}";
  }

Return a hashref with the annotation parsed:

  { chr         => 'field_1',
    start       => 'field_2',
    end         => 'field_3',
    name        => 'field_4',
    score       => 'field_5',
    strand      => 'field_6',
    thick_start => 'field_7',
    thick_end   => 'field_8',
    rgb         => 'field_9'
    blocks      => [ {'size' => 'block size',
                      'start' => 'block start',
                      'end'   => 'block start + block_size',
                      'ref_start' => 'block start on the reference',
                      'ref_end'   => 'block end on the reference'}, ... ],
    seek_pos    => 'Seek position of this line in the file',
  }

=cut

sub bedFileIterator {
  return getFileIterator(file => shift, type => 'bed');
}

=head2 gffFileIterator 

manage GFF3 and GTF2 file format

Example:

  my $it = gffFileIterator($file);
  while (my $annot = $it->()) {
    print "chr    : $annot->{chr}
           start  : $annot->{start}
           end    : $annot->{end}";
  }

Return a hashref with the annotation parsed:

  { chr         => 'field_1',
    source      => 'field_2',
    feature     => 'field_3',
    start       => 'field_4',
    end         => 'field_5',
    score       => 'field_6',
    strand      => 'field_7',
    frame       => 'field_8'
    attributes  => { 'attribute_id' => 'attribute_value', ...},
    seek_pos    => 'Seek position of this line in the file',
  }

gffFileIterator is B<5x faster than Bio-Perl> Bio::Tools::GFF

=cut

sub gffFileIterator {
  my $file = shift;
  my $type = shift;
  return getFileIterator(file => $file, type => $type);
}

=head2 vcfFileIterator 

manage VCF file format

Return a hashref with the annotation parsed:

  { chr => $chr,
    pos     => $pos,
    id      => $id,
    ref     => $ref,
    alt     => [ alt1, alt2, ...],
    qual    => $qual,
    filter  => $filter,
    info    => { AS => value,
                 DP => value,
                 ...
                 ,
  };

=cut

sub vcfFileIterator {
  my $file = shift;
  return getFileIterator(file => $file, type => 'vcf');
}

=head2 chimCTFileIterator

Return a hashref with the chimera parsed:
  
  {
    sample            => $sample,                                                    
    chim_key          => $chim_key,                                                
    chr1              => $chr1,                                                        
    pos1              => $pos1,                                                        
    strand1           => $strand1,                                                  
    chr2              => $chr2,                                                        
    pos2              => $pos2,                                                        
    strand2           => $strand2,                                                  
    chim_value        => $chim_value,                                            
    spanning_junction => $spanning_junction,                              
    spanning_PE       => $spanning_PE,                                          
    class             => $class,                                                      
    comments          => { coment_id => 'comment_value', ... },
    extend_fields     => { extend_field_id => 'extend_field_value', ... },
  }

=cut

sub chimCTFileIterator {
  return getFileIterator(file => shift, type => 'chimCT');
}

=head2 bamFileIterator

BE AWARE this method is only availble if C<samtools> binary is availble.

Return an iterator over a BAM file using a C<samtools view> pipe.

A region can be passed in parameter to restrict the results. In this case
the BAM file must be indexed

Example:

  my $fh = bamFileIterator("file.bam","17:43,971,748-44,105,700");
  while(my $line = <$fh>) {
    my $parsed_line = CracTools::SAMReader::SAMline->new($line);
    // do some stuff
  }

SEE ALSO CracTools::SAMReader::SAMline if you need to parse SAMlines easily

=cut

sub bamFileIterator {
  my ($file,$region) = @_;

  my $fh;
  if($file =~ /\.bam$/) {
    open($fh, "-|", "samtools view $file $region" )or die "Cannot open $file, check if samtools are installed.";
  } else {
    $fh = getReadingFileHandle($file);
  }

  return sub {
    return <$fh>; 
  } 

}

=head2 getSeqFromIndexedRef

=cut

sub getSeqFromIndexedRef {
  my ($ref_file,$chr,$pos,$length) = @_;
  my $region = "$chr:$pos-".($pos+$length);
  return `samtools faidx $ref_file "$region"`;
}

=head1 PARSING LINES

=cut

sub parseBedLine {
  my $line = shift;
  my %args = @_;
  my($chr,$start,$end,$name,$score,$strand,$thick_start,$thick_end,$rgb,$block_count,$block_size,$block_starts) = split("\t",$line);
  my @blocks;

  # Manage blocks if we have a bed12
  if(defined $block_starts) {
    my @block_size = split(",",$block_size);
    my @block_starts = split(",",$block_starts);
    my $cumulated_block_size = 0;
    for(my $i = 0; $i < $block_count; $i++) {
      push(@blocks,{size        => $block_size[$i], 
                   start        => $block_starts[$i], 
                   end          => $block_starts[$i] + $block_size[$i],
                   block_start  => $cumulated_block_size,
                   block_end    => $cumulated_block_size + $block_size[$i],
                   ref_start    => $block_starts[$i] + $start,
                   ref_end      => $block_starts[$i] + $start + $block_size[$i],
                 });
      $cumulated_block_size += $block_size[$i];
    }
  }

  return { chr        => $chr,
    start       => $start, 
    end         => $end, 
    name        => $name,
    score       => $score, 
    strand      => $strand,
    thick_start => $thick_start,
    thick_end   => $thick_end,
    rgb         => $rgb,
    blocks      => \@blocks,
  };
}

=head2 parseGFFLine

=cut

sub parseGFFLine {
  my $line = shift;
  my $type = shift;
  my $attribute_split;
  if($type =~ /gff3/i) {
    $attribute_split = sub {my $attr = shift; return $attr =~ /(\S+)=(.*)/;};
  } elsif ($type eq 'gtf' || $type eq 'gff2') {
    $attribute_split = sub {my $attr = shift; return $attr  =~ /(\S+)\s+"(.*)"/;};
  } else {
    croak "Undefined GFF format (must be either gff3,gtf, or gff2)";
  }
  my($chr,$source,$feature,$start,$end,$score,$strand,$frame,$attributes) = split("\t",$line);
  my @attributes_tab = split(";",$attributes);
  my %attributes_hash;
  foreach my $attr (@attributes_tab) {
    my ($k,$v) = $attribute_split->($attr);
    $attributes_hash{$k} = $v;
  }
  return { chr        => $chr,
    source     => $source,
    feature    => $feature, 
    start      => $start, 
    end        => $end, 
    score      => $score, 
    strand     => $strand,
    frame      => $frame,
    attributes => \%attributes_hash,
  };
}

=head2 parseVCFLine

=cut

sub parseVCFLine {
  my $line = shift;
  my($chr,$pos,$id,$ref,$alt,$qual,$filter,$info) = split("\t",$line);
  my @alts = split(",",$alt);
  my %infos;
  my @infos = split(";",$info);
  foreach (@infos) {
    my($k,$v) = split("=",$_);
    $infos{$k} = $v;
  }

  return { chr => $chr,
    pos     => $pos,
    id      => $id,
    ref     => $ref,
    alt     => \@alts,
    qual    => $qual,
    filter  => $filter,
    info    => \%infos,
  };
}

=head2 parseChimCTLine

=cut

sub parseChimCTLine {                                                     
  my $line = shift;                                                       
  my($id,$name,$chr1,$pos1,$strand1,$chr2,$pos2,$strand2,$chim_value,$spanning_junction,$spanning_PE,$class,$comments,@others) = split("\t",$line);
  my($sample,$chim_key) = split(":",$id);                                 
  my @comments = split(",",$comments);                                    
  my %comments;                                                           
  foreach my $com (@comments){                                            
      my ($key,$val) = split("=",$com);                                   
      $comments{$key} = $val;                                             
  }                                                                       
  my %extend_fields;                                                      
  foreach my $field (@others) {                                           
    my ($key,$val) = split("=",$field);                                   
    if(defined $key && defined $val) {                                    
      $extend_fields{$key} = $val;                                        
    }                                                                     
  }                                                                       
  return {                                                                
    sample            => $sample,                                                    
    chim_key          => $chim_key,                                                
    chr1              => $chr1,                                                        
    pos1              => $pos1,                                                        
    strand1           => $strand1,                                                  
    chr2              => $chr2,                                                        
    pos2              => $pos2,                                                        
    strand2           => $strand2,                                                  
    chim_value        => $chim_value,                                            
    spanning_junction => $spanning_junction,                              
    spanning_PE       => $spanning_PE,                                          
    class             => $class,                                                      
    comments          => \%comments,                                               
    extend_fields     => \%extend_fields,                                     
  };                                                                      
} 


=head1 FILES IO

=head2 getFileIterator

=cut

sub getFileIterator {
  my %args = @_;
  my $file = $args{file};
  my $type = $args{type};
  my $skip = $args{skip};
  my $header_regex = $args{header_regex};
  my $parsing_method = $args{parsing_method};
  my @parsing_arguments = ();


  croak "Missing arguments in getFileIterator" if !defined $file;

  if(!defined $parsing_method && defined $type) {
    if($type =~ /gff3/i || $type eq 'gtf' || $type eq 'gff2') {
      $header_regex = '^#';
      $parsing_method = \&parseGFFLine;
      push (@parsing_arguments,$type);
    } elsif($type =~ /bed/i) {
      $header_regex = '^track\s';
      $parsing_method = \&parseBedLine;
    } elsif ($type =~ /vcf/i) {
      $header_regex = '^#';
      $parsing_method = \&parseVCFLine;
    } elsif ($type =~ /chimCT/i) {
      $header_regex = '^#';
      $parsing_method = \&parseChimCTLine;
    } else {
      croak "Undefined format type";
    }
  }

  # set defaults
  #$header_regex = '^#' if !defined $header_regex;
  $parsing_method = sub{ return { line => shift } } if !defined $parsing_method;

  # We get a filehandle on the file
  my $fh = getReadingFileHandle($file);

  # $curr_pos will hold the SEEK_POS of the current_line
  my $curr_pos = tell($fh);
  # $line will hold the content of the current_line
  my $line = <$fh>;


  # if we need to skip lines
  if(defined $skip) {
    for(my $i = 0; $i < $skip; $i++) {
      $curr_pos = tell($fh);
      $line = <$fh>;
    }
  }
  
  # Skip line that match a specific regex
  if(defined $header_regex) {
    while($line =~ /$header_regex/) {
      $curr_pos = tell($fh);
      $line = <$fh>;
    }
  }

  # The iterator itself (an unnamed subroutine
  return sub {
    if (defined $line) {
      chomp $line;
      # We parse the line with the appropriate methd
      my $parsed_line = $parsing_method->($line,@parsing_arguments);
      # Add seek pos to the output hashref
      $parsed_line->{seek_pos} = $curr_pos;
      $curr_pos = tell($fh);
      $line = <$fh>; # Get next line
      return $parsed_line;
    } else {
      return undef;
    }
  };
}

=head2 getReadingFileHandle

Return a file handle for the file in argument.
Display errors if file cannot be oppenned and manage gzipped files (based on .gz file extension)

Example:
  
  my $fh = getReadingFileHandle('file.txt.gz');
  while(<$fh>) {
    print $_;
  }
  close $fh;

=cut

sub getReadingFileHandle {
  my $file = shift;
  my $fh;
  if($file =~ /\.gz$/) {
    open($fh,"gunzip -c $file |") or die ("Cannot open $file");
  } else {
    open($fh,"< $file") or die ("Cannot open $file");
  }
  return $fh;
}

=head2 getWritingFileHandle

Return a file handle for the file in argument.
Display errors if file cannot be oppenned and manage gzipped files (based on .gz file extension)

Example:
  
  my $fh = getWritingFileHandle('file.txt.gz');
  print $fh "Hello world\n";
  close $fh;

=cut

sub getWritingFileHandle {
  my $file = shift;
  my $fh;
  if($file =~ /\.gz$/) {
    open($fh,"| gzip > $file") or die ("Cannot open $file");
  } else {
    open($fh,"> $file") or die ("Cannot open $file");
  }
  return $fh;
}

=head2 getLineFromSeekPos 

  getLineFromSeekPos($filehandle,$seek_pos);

return a chomped line at a seeking position.

=cut

sub getLineFromSeekPos {
  my($fh,$seek_pos) = @_;
  seek($fh,$seek_pos,SEEK_SET);
  my $line = <$fh>;
  chomp($line);
  return $line;
}

1;
