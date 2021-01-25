#!/usr/bin/perl -w

#########################
#
#    Bubblyhelix - reversify.pl
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

# usage: perl reversify.pl <fasta_directory>

use File::Path 'mkpath';

my $dir=$ARGV[0];
die unless (-d $dir); # input arg must be a directory

mkpath(["$dir/revcompl"],0,0777); # make subdir called revcompl (unless it exists)

my %compl=(
# taken from http://www.chem.qmul.ac.uk/iubmb/misc/naseq.html
  A => "T",
  T => "A",
  C => "G",
  G => "C",
  N => "N",
  S => "S",
  W => "W",
  B => "V",
  V => "B",
  D => "H",
  H => "D",
  K => "M",
  M => "K",
  R => "Y",
  Y => "R",
# preserves upper/lower case:
  a => "t",
  t => "a",
  c => "g",
  g => "c",
  n => "n",
  s => "s",
  w => "w",
  b => "v",
  v => "b",
  d => "h",
  h => "d",
  k => "m",
  m => "k",
  r => "y",
  y => "r",
);

my @fa_files=glob "$dir/*.fa"; # all files ending with .fa

while (@fa_files) {
  my $forwardfile=shift @fa_files;
  my @path=split "/", $forwardfile;
  my $file=pop @path;
  my $revfile=join "/",(@path,"revcompl",$file);
  next if (-e $revfile); # skip those that already are "reversified"
  
  print "creating $revfile\n";
  
  # read fastafile:
  open(FA, $forwardfile) ||
    die("Cannot open $forwardfile: $!");
  chomp(my $header=<FA>);
  chomp(my @lines=<FA>);
  close FA;
  
  # write reverse fastafile:
  open RC, ">$revfile"
    or die "cannot create $revfile: $!";
  print RC "> the reverse complement of [ $header ]\n";
  while (@lines) {
    # string converted to array and back to string:
    my @revstring;
    foreach $char (split "", pop @lines) {
      if (exists $compl{$char}) {
        unshift @revstring, $compl{$char};
      } else {
        warn "$char unknown";
      }
    }
    print RC join "", @revstring;
  }  
  close RC;
}

