#!/usr/bin/perl -w

#########################
#
#    Bubblyhelix - random_fasta.pl
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

# usage: perl random_fasta.pl <length> > <filename.fa>
#
# This script produces a fastafile with a random sequence of given length,
# but without the structure of biological sequences (exons, introns, isochores, etc.)
# that can induce structure in the melting profiles. 


$N=$ARGV[0];

print ">\n"; # fasta header
$x=0;
@base=qw/A C G T/;
while ($x<$N) {
  foreach (1..80) {
    print $base[int(rand(4))];
    $x++;
    last if ($x==$N);
  }
  print "\n";    
}  