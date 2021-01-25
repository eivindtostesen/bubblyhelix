#########################
#
#    Bubblyhelix - Melting_strand_region.pm
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

package Melting_strand_region;

@ISA=qw(Melting_strand);
use Melting_strand;

sub Z01 { $_[0]->{"Z01"}->[$_[1] - $_[0]->{"pf array offset"}] }
sub Z11 { $_[0]->{"Z11"}->[$_[1] - $_[0]->{"pf array offset"}] }
sub Z10 { $_[0]->{"Z10"}->[$_[1] - $_[0]->{"pf array offset"}] }

sub xpfvector_to_arrays {
  # inserts xpfvector in arrays at position x
  my $strand=shift;
	my ($x,$Z01,$Z11,$Z10,$Zeo)=@_[0,4,5,6,2]; # @_ is xpfvector
  $strand->{"Z01"}->[$x - $strand->{"pf array offset"}]=$Z01;
  $strand->{"Z11"}->[$x - $strand->{"pf array offset"}]=$Z11;
  $strand->{"Z10"}->[$x+1 - $strand->{"pf array offset"}]=$Z10;
  $strand->{"Zeo values"}->[$x/$strand->sw->modulus_Z]=$Zeo;
}
sub forward_and_reverse {
  my $strand=shift;
  my ($start,$end)=@_;

  # fill arrays both ways in the region:
  my @xpfvector1=$strand->pf_fill_arrays($start,$end);
  my @xpfvector2=$strand->rev->pf_fill_arrays( $strand->revloc($start,$end) ); 

  # compute Ztotal:
  if ($start==1 and $end==$strand->length) {
    # then both strands get Ztotal from their own xpfvector
    $strand->set_Ztotal($strand->compute_Ztotal(@xpfvector1));
    $strand->rev->set_Ztotal($strand->rev->compute_Ztotal(@xpfvector2));
  } elsif ($start-1>=$strand->length-$end) {
    # then both strands get Ztotal from the forward xpfvector1
    if ($end<$strand->length) {
      # post sweep from $end to strand length:
	    @xpfvector1=$strand->pf_recursion_no_arrays($strand->length, @xpfvector1);
    }  
    $strand->set_Ztotal($strand->compute_Ztotal(@xpfvector1));
    # copy Ztotal to reverse strand:
    $strand->rev->set_Ztotal($strand->Ztotal);
  } else {
    # then both strands get Ztotal from the reverse xpfvector2
    if ($start>1) {
      # post sweep from $end to strand length:
	    @xpfvector2=$strand->rev->pf_recursion_no_arrays($strand->length, @xpfvector2);
    }  
    $strand->rev->set_Ztotal($strand->rev->compute_Ztotal(@xpfvector2));
    # copy Ztotal to forward strand:
    $strand->set_Ztotal($strand->rev->Ztotal);
  }  
}
sub pf_fill_arrays {
  my $strand=shift;
  my ($start,$end)=@_;

  # erase previous arrays:
  $strand->tabula_rasa;

  # initialize:
  $strand->{"pf array offset"}=$start;
  $strand->sw->set_iterator(1); # rewind
	my @xpfvector=$strand->pf_first_xpfvector;

  if ($start>1) {
    # pre sweep from 1 to $start:
    @xpfvector=$strand->pf_recursion_no_arrays($start, @xpfvector);
  }  
  
  # fill arrays:
  $strand->xpfvector_to_arrays(@xpfvector); # here x=$start in @xpfvector
	@xpfvector=$strand->pf_recursion(
    sub {
      push (@{ $strand->{"Z01"} }, $_[0]);
      push (@{ $strand->{"Z11"} }, $_[1]);
      push (@{ $strand->{"Z10"} }, $_[2]);
    },
    sub {
      push (@{ $strand->{"Zeo values"} }, $_[0]); # push on the staircase
    },
    $end,
    @xpfvector,
  );
  
	$strand->calculate_Zeo_jumps($start,$end);
  return @xpfvector;
}
sub pf_recursion_no_arrays {
  my $strand=shift;
  $strand->pf_recursion(sub {}, sub {}, @_);
} 

1;


