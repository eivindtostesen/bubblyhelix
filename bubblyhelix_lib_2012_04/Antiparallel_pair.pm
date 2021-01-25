#########################
#
#    Bubblyhelix - Antiparallel_pair.pm
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

package Antiparallel_pair;

sub new_pair {
  my $class=shift @_;
  my $strand1=$class->new;
  my $strand2=$class->new;
  $strand1->pair_with($strand2);
  $strand1->set_pair_length(shift) if (@_); # sets both!
  return ($strand1,$strand2);
}
sub rev { $_[0]->{"partner"} } # reverse complement
sub pair_with {
  my ($strand1,$strand2)=@_;
  $strand1->unpair;
  $strand2->unpair;
  $strand1->{"partner"}=$strand2;
  $strand2->{"partner"}=$strand1;
}
sub is_paired {
  my $strand=shift;
  exists $strand->{"partner"};
}
sub unpair {
  my $strand=shift;
  if ($strand->is_paired) {
    delete $strand->rev->{"partner"};
    delete $strand->{"partner"};
  }
}

sub length { $_[0]->{"length"} }
sub set_pair_length {
  my $strand=shift;
  $strand->{"length"}=$_[0];
  $strand->rev->{"length"}=$_[0] if ($strand->is_paired);
}
sub revpos {	$_[0]->length+1-$_[1] }
sub revloc {
  my $strand=shift;
	if (@_ == 1) {
	  # this means the NN at [ $_[0]-1, $_[0] ]
		return $strand->revpos($_[0]) + 1;
	} elsif (@_ > 1) {
    return map {$strand->length+1-$_} reverse @_;
	}
}



1;


