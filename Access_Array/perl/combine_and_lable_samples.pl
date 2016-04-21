#!/usr/bin/perl -w

use strict;
use warnings;
use PerlIO::gzip;
use File::Basename;


while (<>){
  my $filename = basename ($ARGV, ".fastq");
  if($. % 4 == 1){
    print  "@"."RECORD:".$filename."_".$.."\n";
  }
  else {
    print $_;
  }
}
