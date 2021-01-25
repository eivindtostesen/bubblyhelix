#!/usr/bin/perl -w

#########################
#
#    Bubblyhelix - pfdbm2txt.pl
#
#    Copyright 2012 Eivind Tostesen.
#
#    This file is part of Bubblyhelix.
#
#    Bubblyhelix is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    Bubblyhelix is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with Bubblyhelix.  If not, see <http://www.gnu.org/licenses/>.
#
#########################

# This script dumps a pfdbm directory (containing DBM files) to a textfile-based "copy".
# Unlike DBM files, textfiles are human-readable and portable between machines.
#
# usage: perl pfdbm2txt.pl <pfdbm_source/> <pfdbm_targetdir> > <pfdbm2txt_log.txt> 2>&1

use File::Copy;
use DB_File;

&copyrecursive($ARGV[0],$ARGV[1]);

sub copyrecursive {
  my ($sourcedir,$targetdir)=@_;
  die "$targetdir already exists!" if (-e $targetdir);
  if (-d $sourcedir) {
    mkdir $targetdir, 0755 or die;
    my $sd;
    opendir $sd, $sourcedir or die;
    foreach $name (readdir($sd)) {
      next if $name eq "." or $name eq "..";
      my ($orig,$copy)=("$sourcedir/$name","$targetdir/$name");
      if (-d $orig) {
        &copyrecursive($orig,$copy);
      } elsif ($name =~ /backward|forward|reverse|Ztotal_b|Ztotal_f|Ztotal_forward|Ztotal_reverse/) {
        &pfdbm_dumper($orig,"$copy.txt");
      } elsif (-f $orig) {
        !system "cp -p $orig $copy" or die; # use this option on Mac OS X
        #!system "cp -a $orig $copy" or die; # use this option on Linux
      } else {
        print "$orig was not copied!\n";
      }
    }
    close $sd;
  } else {
    die "$sourcedir is not a directory!";
  }  
}
sub pfdbm_dumper {
  # dumps a DBM file to a textfile
  my ($dbm,$txt)=@_;
  my $dbmcopy=$txt."_$$.tmp";
  my %STATE;
  # make temporary copy:
  copy($dbm,$dbmcopy) or die "Copy failed: $!";
  # opening the copy keeps the time stamp of the original:
  dbmopen(%STATE, "$dbmcopy", 0644)
    or die "cannot dbmopen $dbmcopy: $!";
  open my $new, ">", $txt or die "cannot create $txt: $!"; 
  foreach $x (sort {$a<=>$b} keys %STATE) {
    print {$new} "$x\t$STATE{$x}\n";
  }
  close $new;
  dbmclose(%STATE);
  unlink $dbmcopy;
}