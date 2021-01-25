#########################
#
#    Bubblyhelix - Iupac.pm
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

# Provides averaging for IUPAC symbols other than ACGT

package Iupac;

my @known=qw/A C G T/;

sub include {
  my ($s_11,$s_1,$s_010,$unknown)=@_;

  # calculate the "unknown" s1 and s010 values:
  foreach $Z (@$unknown) {
    # init:
    $s_1->{$Z}=1;
    $s_010->{$Z}=1;
	# geometric means:
    foreach $X (@known) {
      $s_1->{$Z}*=$s_1->{$X}**(1/4);
      $s_010->{$Z}*=$s_010->{$X}**(1/4);
    }
  }

  # calculate the "known-unknown" and "unknown-known" s11-values:
  foreach $Z (@$unknown) {
    # init:
	foreach $X (@known) {
	  $s_11->{$Z}->{$X}=1;
	  $s_11->{$X}->{$Z}=1;
	}
	# geometric means:
    foreach $X (@known) {
      foreach $Y (@known) {
        $s_11->{$Z}->{$Y}*=$s_11->{$X}->{$Y}**(1/4);
        $s_11->{$X}->{$Z}*=$s_11->{$X}->{$Y}**(1/4);
      }
    }
  }
  
  # calculate the "unknown-unknown" s11-values:
  foreach $Z1 (@$unknown) {
    foreach $Z2 (@$unknown) {
      # init:
	  $s_11->{$Z1}->{$Z2}=1;
      # geometric mean:
      foreach $X (@known) {
	    foreach $Y (@known) {
	      $s_11->{$Z1}->{$Z2}*=$s_11->{$X}->{$Y}**(1/16);
	    }
	  }
    }
  }
}

1;

