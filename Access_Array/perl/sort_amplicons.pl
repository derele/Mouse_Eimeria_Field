#!/usr/bin/perl -w
use strict;
use warnings;
use Data::Dumper;

## development version... might change with time

die "USAGE: sort_amplicons.pl primer_csv Read1_file Read2_file output_folder\n" unless scalar(@ARGV) == 4;

## csv file should have 4 columns: Primer Name F, R, Primer seq F, R
## the two read files have to be paired end fastq format 
## folder will be created if it exists files will be overwritten

my ($primerF, $primerR, $nameF, $nameR) = IUPAC_Translator($ARGV[0]);
my $R1_file = $ARGV[1];
my $R2_file = $ARGV[2];
my $out_folder = $ARGV[3];

## create a hash of arrays for each unique forward primer and the
## reverse primers used in combination with it
my %amps;
for (my $i=0; $i<scalar(@$primerF); $i++ ) {
  push( @{$amps{$$primerF[$i]}},  $$primerR[$i] );
} 

my %FQ = read_fastq_both_reads(); ## the input sequences 
my %seqALL; ## the output data structure
my @match; ## a temporary store for matches

print STDERR scalar( (keys %FQ)) . " fastq records read \n";

## the sorting machine for good matches
foreach  my $id (keys %FQ) {
  foreach my $primerX (keys %amps) { # based on forward 
    my $primerXq=$primerX; # a quoted version for the hash access
    $primerX=~ s/\"//g; # a unquoted version for the regex    
    if ($FQ{$id}{seqF} =~ s/^$primerX//) { # a forward primer is found
      foreach my $primerY ( @{ $amps{$primerXq} }){       
        $primerY=~ s/\"//g; # a unquoted version for the regex
        if ($FQ{$id}{seqR} =~ s/^$primerY//){ # a reverse primer is found
          $seqALL{$primerX}{$primerY}{$id} = $FQ{$id}; 
          ## truncate the quality lines accordingly
          $seqALL{$primerX}{$primerY}{$id}{qualF} = 
            substr($seqALL{$primerX}{$primerY}{$id}{qualF}, 
                   -length($seqALL{$primerX}{$primerY}{$id}{seqF}));
          $seqALL{$primerX}{$primerY}{$id}{qualR} = 
            substr($seqALL{$primerX}{$primerY}{$id}{qualR}, 
                   -length($seqALL{$primerX}{$primerY}{$id}{seqR}));; 
          ## from here: full amps sorting, still way to complicated...
          my $RCrev=revComp($primerY);
          if ($seqALL{$primerX}{$primerY}{$id}{seqF} =~ s/$RCrev.*//){
            my $FCrev=revComp($primerX);
            $seqALL{$primerX}{$primerY}{$id}{seqR} =~ s/$FCrev.*//;
            ## truncate quality lines accordingly
            $seqALL{$primerX}{$primerY}{$id}{qualF} = 
              substr($seqALL{$primerX}{$primerY}{$id}{qualF},0,
                     length($seqALL{$primerX}{$primerY}{$id}{seqF}));
            $seqALL{$primerX}{$primerY}{$id}{qualR} = 
              substr($seqALL{$primerX}{$primerY}{$id}{qualR},0,
                     length($seqALL{$primerX}{$primerY}{$id}{seqR}));
            ## add a special identifier sequences
            $seqALL{$primerX.":full"}{$primerY.":full"}{$id} = 
              $seqALL{$primerX}{$primerY}{$id};
            ## delete the original
            delete ($seqALL{$primerX}{$primerY}{$id});
          }
          push @match, $id; 
        }
      }
    }
  }
}


## remove all sequences for which we had matches from the original
## fastq hash
delete @FQ{@match};
print STDERR scalar( (keys %FQ)) . " fastq records not sorting into amplicons\n";

## the forward matching but NOT reverse matching primers
foreach my $id (keys %FQ) { # %FQ only reads that failed before     
  foreach my $primerX (keys %amps) { # based on forward 
    my $primerXq=$primerX; # a quoted version for the hash access
    $primerX=~ s/\"//g; # a unquoted version for the regex    
    if ($FQ{$id}{seqF} =~ s/^$primerX//) { # a forward primer is found
      $seqALL{$primerX}{failR}{$id} = $FQ{$id}; 
      push @match, $id; 
    }
  }
}

## remove again all sequences for which we had matches from the original
## fastq hash
delete @FQ{@match};

## NOT forward but reverse primer matching  
my @primerRu = uniq(@$primerR);
foreach  my $id (keys %FQ) {  
  foreach my $primerY (@primerRu) {
    my $primerYq=$primerY;
    $primerY=~ s/\"//g; # a unquoted version for the regex    
    if ($FQ{$id}{seqR} =~ s/^$primerY//) { # a reverse primer is found
      $seqALL{failF}{$primerYq}{$id} = $FQ{$id};
      push @match, $id; 
    }
  }
}

delete @FQ{@match};

## left are only the neither forward nor reverse matching
foreach  my $id (keys %FQ) { 
  $seqALL{failF}{failR}{$id} = $FQ{$id};
}

## traverse the sorted structures and print 
foreach my $F (keys %seqALL) {
  foreach my $R (keys %{$seqALL{$F}} ) {
    my $c = (keys %{$seqALL{$F}{$R}} );
    my $pNameF = primerSeq2Name($F);
    $pNameF =~ s/\"//g;
    my $pNameR = primerSeq2Name($R);
    $pNameR =~ s/\"//g;
    ## write a simple tsv table
    print ($pNameF."\t".$pNameR."\t".$F."\t".$R."\t".$c."\n");
    ## write two fastq sequences for each amplicon
    write_fastq_both_reads($seqALL{$F}{$R},
                           $pNameF."_".$pNameR,
                           $F."_".$R);
  }
}


############################# subroutines ###########################

## read fastq from two files
sub read_fastq_both_reads {
  open R1, "<", $R1_file or die; 
  open R2, "<", $R2_file or die; 
  my %FQF;
  chomp(my @fastqF = <R1>);
  chomp(my @fastqR = <R2>);
  for  (my $i = 0; $i< scalar (@fastqF); $i=$i+4) {
    $FQF{$fastqF[$i]}{"headerR"} = $fastqR[$i];
    $FQF{$fastqF[$i]}{"seqF"} = $fastqF[$i+1];
    $FQF{$fastqF[$i]}{"seqR"} = $fastqR[$i+1];
    $FQF{$fastqF[$i]}{"qualF"} = $fastqF[$i+3];
    $FQF{$fastqF[$i]}{"qualR"} = $fastqR[$i+3];
  }
  close R1; 
  close R2;
  return(%FQF)
}

## write fastq to two files
sub write_fastq_both_reads {
  my ($R, $amplicon, $primcat) = (@_);
  my %R = %$R; ## We expect a hash ref here
  unless (-e $out_folder or mkdir $out_folder) {
    die "Unable to create $out_folder: $!"};
  open R1, ">", $out_folder."/".
    $amplicon."_R1.fastq" or die; 
  open R2, ">", $out_folder."/".
    $amplicon."_R2.fastq" or die; 
  foreach my $id (keys %R) {
    print R1 $id.":".$amplicon."\n".
      $R{$id}{seqF}."\n".
        "+PRIMER:".$primcat."\n".
          $R{$id}{qualF}."\n";
    print R2 $R{$id}{headerR}.":".$amplicon."\n".
      $R{$id}{seqR}."\n".
        "+PRIMER:".$primcat."\n".
          $R{$id}{qualR}."\n";
  }
}

## Tranlate primer from csv to regex 
sub IUPAC_Translator {
  my $input_primer_csv = shift; # argument: file path to csv
  open my $primer_csv, '<', $input_primer_csv
    or die "Can't read from $input_primer_csv: $!";
  my %symbols_to_sub = (
                        "V" => "(A|C|G)",
                        "D" => "(A|T|G)",
                        "B" => "(T|G|C)",
                        "H" => "(A|T|C)",
                        "W" => "(A|T)",
                        "S" => "(C|G)",
                        "K" => "(T|G)",
                        "M" => "(A|C)",
                        "Y" => "(C|T)",
                        "R" => "(A|G)",
                        "N" => "(A|T|C|G)",
                       );
  
  ## read the primers
  my (@nameF, @nameR, @primerF, @primerR );
  while (<$primer_csv>) {
    chomp();
    my @split = split(';');
    next if $. < 2;             #skip header line
    push @nameF, $split[0];
    push @nameR, $split[1];
    push @primerF, $split[2];
    push @primerR, $split[3];
  }

  ## do the translation int regular expressions
  foreach my $primer ( @primerF, @primerR ) {
    foreach my $key ( keys %symbols_to_sub ) {
      $primer =~ s/$key/$symbols_to_sub{$key}/g;
    }
  }
  return( \@primerF, \@primerR, \@nameF, \@nameR)
}

sub uniq {
  my %seen;
  grep !$seen{$_}++, @_;
}


## reverse complement (for regex)
sub revComp{
    my $seq = shift(@_);
    $seq    = reverse($seq);
    $seq   =~ tr/A)TC(G/T(AG)C/;
    return($seq);
}


sub primerSeq2Name{
  my $S = shift(@_);            # just a sequence
  for (my $i=0; $i<scalar(@$nameF); $i++ ) {
    my $PrF = $$primerF[$i];    # using global vars for ...
    $PrF =~ s/\"//g;
    if ( $PrF eq $S ) {
      return ($$nameF[$i]);     # everything else 
    } 
    elsif ($PrF.":full" eq $S) {
      return ($$nameF[$i].":full");
    }
  }
  for (my $i=0; $i<scalar(@$nameR); $i++ ) {
    my $PrR = $$primerR[$i];
    $PrR =~ s/\"//g;
    if ($PrR eq $S) {
      return ($$nameR[$i]);
    } 
    elsif ($PrR.":full" eq $S) {
      return ($$nameR[$i].":full");
    }
  }
  return ("fail");              # if we get to here we have no primer
}
