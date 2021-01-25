#########################
#
#    Bubblyhelix - Melting_strand.pm
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

package Melting_strand;

@ISA=qw(Antiparallel_pair);
use Antiparallel_pair;

use POSIX (); # empty import
sub round_et {
  return ($_[0] >= 0) ? POSIX::floor( $_[0]+0.5) : POSIX::ceil( $_[0]-0.5);
}

sub new {
  my $class=shift @_;
  bless { @_ }, $class;
};
sub empty_cache { $_[0]->{"cache"}=() }
sub sw { $_[0]->{"seq 2 sw"} }
sub set_sw { $_[0]->{"seq 2 sw"}=$_[1] }
sub Z01 { $_[0]->{"Z01"}->[$_[1]] }
sub Z11 { $_[0]->{"Z11"}->[$_[1]] }
sub Z10 { $_[0]->{"Z10"}->[$_[1]] }
### Total
sub Ztotal { @{ $_[0]->{"Ztotal"} } }
sub divide_by_Ztotal {	($_[1]/$_[0]->{"Ztotal"}->[0], $_[2]-$_[0]->{"Ztotal"}->[1]) }
sub set_Ztotal {
  my $strand=shift;
	die unless (@_);
  $strand->{"Ztotal"}=[@_]; # e.g. from database or from strand->rev
}
sub compute_Ztotal {
  my $strand=shift;
   # @_ is xpfvector:
  die unless ($_[0]==$strand->length);
  ($strand->rev->sw->sigma_tail * ($_[3]-$_[6]) + $_[6] , $_[2]);
}
### eo
sub Zeo { $_[0]->{"Zeo values"}->[$_[1]/$_[0]->sw->modulus_Z] } # Piecewise Constant Exponent Offset (PCEO)
sub calculate_Zeo_jumps {
  my $strand=shift;
  my ($start,$end)=@_;
  $strand->{"Zeo jumps"}={ map {$_*$strand->sw->modulus_Z,$strand->{"Zeo values"}->[$_]-$strand->{"Zeo values"}->[$_-1]} 1+int($start/$strand->sw->modulus_Z)..int($end/$strand->sw->modulus_Z) };
}
sub sum_pceo {
  # PCEO scheme for summing f(x) (with Zeo) over x in [start,end]
  my $strand=shift;
  my ($start,$end,$f)=@_; # $f is a sub ref
  my $sum=0;
  foreach $x ($start..$end) {
	  $sum/=10**$strand->{"Zeo jumps"}->{$x} if ($strand->{"Zeo jumps"}->{$x});
	  $sum+=$f->($x);
  }
  ($sum, $strand->Zeo($end));  # result = $sum*10**( $strand->Zeo($end) )
}
sub sum_eo {
  # EO scheme for summing (unordered) stream of (fp,eo)-pairs
  my $strand=shift;
  my ($iterator)=@_;
	my (%sums,@i,$eo,$tot);
	$sums{$i[1]}+=$i[0] while @i=&$iterator;
	foreach (sort {$a<=>$b} keys %sums) {
	  $tot/=10**($_-$eo) if (defined $eo); # skip the first time
		$tot+=$sums{$_};
		$eo=$_;
	}
	return ($tot,$eo);
}
####################### PF

sub tabula_rasa {
  my $strand=shift;
  $strand->empty_cache;
  $strand->{"Z01"}=[];
  $strand->{"Z11"}=[];
  $strand->{"Z10"}=[];
  $strand->{"Zeo values"}=[];
}
sub xpfvector_to_arrays {
  # inserts xpfvector in arrays at position x
  my $strand=shift;
	my ($x,$Z01,$Z11,$Z10,$Zeo)=@_[0,4,5,6,2]; # @_ is xpfvector
  $strand->{"Z01"}->[$x]=$Z01;
  $strand->{"Z11"}->[$x]=$Z11;
  $strand->{"Z10"}->[$x+1]=$Z10;
  $strand->{"Zeo values"}->[$x/$strand->sw->modulus_Z]=$Zeo;
}
sub forward_and_reverse {
  my $strand=shift;
  
  # fill arrays both ways
  my @xpfvector1=$strand->pf_fill_arrays;
  my @xpfvector2=$strand->rev->pf_fill_arrays; 
  
  # compute Ztotal:
  $strand->set_Ztotal($strand->compute_Ztotal(@xpfvector1));
  $strand->rev->set_Ztotal($strand->rev->compute_Ztotal(@xpfvector2));
}
sub pf_fill_arrays {
  my $strand=shift;

  # erase previous arrays:
  $strand->tabula_rasa;

  # initialize:
  $strand->sw->set_iterator(1); # rewind
	my @xpfvector=$strand->pf_first_xpfvector;

  # fill arrays:
  $strand->xpfvector_to_arrays(@xpfvector); # here x=1 in @xpfvector
	@xpfvector=$strand->pf_recursion(
    sub {
      push (@{ $strand->{"Z01"} }, $_[0]);
      push (@{ $strand->{"Z11"} }, $_[1]);
      push (@{ $strand->{"Z10"} }, $_[2]);
    },
    sub {
      push (@{ $strand->{"Zeo values"} }, $_[0]); # push on the staircase
    },
    $strand->length,
    @xpfvector,
  );

	$strand->calculate_Zeo_jumps(1,$strand->length);
  return @xpfvector;
}
sub pf_first_xpfvector {
  # returns xpfvector at x=1
  my $strand=shift;
  return (
    1, # 0 start
    $strand->sw->modulus_Z, # 1 next_carry
    0, # 2 Zeo
    $strand->sw->s010(1), # 3 partial_sum_Z10
    1, # 4 Z01
    0, # 5 Z11
    $strand->sw->s010(1), # 6 Z10
    $strand->sw->sigma_tail, # 7 L_tail
    (map 0, 0..$strand->sw->C_last), # 8,9,... W
  );
}
sub pf_recursion {
  # computes partition functions by bottom-up dynamic programming (recursion)
  my $strand=shift;
  my ($callback1,$callback2,$end,$start,$next_carry,$Zeo,$partial_sum_Z10,$Z01,$Z11,$Z10,$L_tail,@W)=@_;
  foreach $x ($start+1..$end) { # the main loop!
    my ($s011,$s110,$s11,$s010)=$strand->sw->next_transfer_matrix;  # positions ($x-1,$x)
    $Z11=$s11*($s011*$Z01+$Z11);
    $Z01=$L_tail;
    foreach $m (0..$#W) {
      $Z01+=$W[$m];
    }
    foreach $m (0..$#W) {
      $W[$m]=$Z10*$strand->sw->C1($m)+$W[$m]*$strand->sw->C2($m);
    }
    $Z10=$s110*$Z11+$s010*$Z01; # = $strand->Z10($x+1)
    $partial_sum_Z10+=$Z10;
    if ($x==$next_carry) {  #carry at the modulus
      my $carry=&round_et( log($partial_sum_Z10)/log(10));
      foreach ($partial_sum_Z10, $Z01, $Z11, $Z10, $L_tail, @W) {
        $_/=10**$carry;
      }
      $Zeo+=$carry;
      $next_carry=&round_et( $next_carry+$strand->sw->modulus_Z ); # prepare next integer carry_position
      $callback2->($Zeo);
    }
    $callback1->($Z01,$Z11,$Z10);
  }
  return ($end,$next_carry,$Zeo,$partial_sum_Z10,$Z01,$Z11,$Z10,$L_tail,@W);
}
################# PROBABILITIES
sub p_1 {
  my $strand=shift;
  my $x=$_[0];
  if ($x>1 and $x<$strand->length) {
    my $rx=$strand->revpos($x);
    my $p=( $strand->Z01($x) * $strand->sw->s1($x) * $strand->rev->Z11($rx)
					 +$strand->Z11($x) * $strand->sw->s1($x) * $strand->rev->Z01($rx)
					 +$strand->Z01($x) * $strand->sw->s010($x) * $strand->rev->Z01($rx)
					 +$strand->Z11($x) * $strand->rev->Z11($rx) );
    my $eo=$strand->Zeo($x)+$strand->rev->Zeo($rx);
		return $strand->divide_by_Ztotal($p,$eo);
  } elsif ($x==1) {
    return $strand->p_left(1);
  } elsif ($x==$strand->length) {
    return $strand->rev->p_left(1);
  }  
}
sub p_right {
  my $strand=shift;
  my $x=$_[0];
  if ($x>=1 and $x<$strand->length) {
  	return $strand->divide_by_Ztotal($strand->Z10($x+1)*$strand->rev->sw->sigma_tail, $strand->Zeo($x));
  } elsif ($x==$strand->length) {
  	return $strand->divide_by_Ztotal($strand->Z10($x+1), $strand->Zeo($x));
	}
}
sub p_01 {
  my $strand=shift;
  my $x=$_[0];
  if ($x>=2 and $x<=$strand->length) {
    my $rx=$strand->revpos($x);
    my $p=$strand->Z01($x)*$strand->rev->Z10($rx+1);
    my $eo=$strand->Zeo($x)+$strand->rev->Zeo($rx);
		return $strand->divide_by_Ztotal($p,$eo);
  }
}
sub p_left { $_[0]->rev->p_right( $_[0]->revpos($_[1]) ) }
sub p_10 { $_[0]->rev->p_01( $_[0]->revloc($_[1]) ) } # revloc of nn (not revpos)
sub p_bubble {
  my $strand=shift;
  my ($x,$y)=@_;
  if ($x+1<$y and $x>=1 and $y<=$strand->length) {
    my $ry=$strand->revpos($y);
    my $p=$strand->Z10($x+1)*$strand->sw->omega($y-$x)*$strand->rev->Z10($ry+1);
    my $eo=$strand->Zeo($x)+$strand->rev->Zeo($ry);
		return $strand->divide_by_Ztotal($p,$eo);
  }
}
############### PROBABILITY PROFILES
sub p_1_profile {
  my $strand=shift;
  my ($start,$end,$output)=@_;
  foreach $i ($start..$end) {
    $output->($strand->p_1($i));
  }
}
################ PROBABILITY SUMS AND BOUNDS
sub helicity {
  my $strand=shift;
  my ($start,$end)=(1,$strand->length);
  ($start,$end)=@_ if (@_);
	my $i=$start-1; # rewind
  my ($sump1,$sump1_eo)=$strand->sum_eo( sub { ++$i<=$end  ?  $strand->p_1($i)  :  () } );
	($sump1/($end-$start+1), $sump1_eo);
}
sub sum_loc_Z10_cache {
  #returns (from cache or calculation) the sum of Z10 over input interval [start,end] (+1)
  my $strand=shift;
	unless (exists $strand->{"cache"}->{"sum loc Z10"}->{"@_"}) {
    $strand->{"cache"}->{"sum loc Z10"}->{"@_"}=[ $strand->sum_pceo(@_, sub { $strand->Z10($_[0]+1) } ) ];
  }
  @{ $strand->{"cache"}->{"sum loc Z10"}->{"@_"} };    
}
sub sum_p_right {
  # calculates the sum of p_right over input interval [start,end]
  my $strand=shift;
  my ($start,$end)=@_;
	if ($end<$strand->length) {
    my ($sumZ10,$sumZ10_eo)=$strand->sum_loc_Z10_cache($start,$end);
	  return $strand->divide_by_Ztotal( $sumZ10 * $strand->rev->sw->sigma_tail, $sumZ10_eo);
	} elsif ($end==$strand->length) {
  	my $i=$start-1; # rewind
	  return $strand->sum_eo( sub { ++$i<=$end  ?  $strand->p_right($i)  :  () } );
  }
}
sub sum_p_left { $_[0]->rev->sum_p_right( $_[0]->revloc(@_[1,2]) ) }
sub sum_p_bubble {
  # calculates the sum of p_bubble over input frame [x1,x2]x[y1,y2]
  # with speed adaptation
  my $strand=shift;
  my ($x1,$x2,$y1,$y2)=@_;
  if (($x2-$x1)*($y2-$y1)>($strand->sw->C_last+1)*($x2-$x1+$y2-$y1)) { #if big frame
    $strand->sum_p_bubble_fast($x1,$x2,$y1,$y2); # multiexp is fastest
  } else { #if small frame
    $strand->sum_p_bubble_slow($x1,$x2,$y1,$y2); # full summation is fastest
  }
}
sub sum_p_bubble_slow {
  # calculates the sum of p_bubble over input frame [x1,x2]x[y1,y2]
  my $strand=shift;
  my ($x1,$x2,$y1,$y2)=@_;
  my @ry=$strand->revloc($y1,$y2);
  my ($inner_sum,$inner_sum_eo);
  my ($outer_sum,$outer_sum_eo)=$strand->sum_pceo(
    $x1,$x2,
	  sub {
	    my $x=$_[0];
	    my $rx=$strand->revpos($x);
	    ($inner_sum,$inner_sum_eo)=$strand->rev->sum_pceo(
			  @ry,
		    sub { $strand->sw->omega($rx-$_[0]) * $strand->rev->Z10($_[0]+1) },
	    );
	    return $strand->Z10($x+1) * $inner_sum;
	  },
  );
  $strand->divide_by_Ztotal($outer_sum, $inner_sum_eo+$outer_sum_eo);
}
sub sum_p_bubble_fast {
  # calculates the sum of p_bubble over input frame [x1,x2]x[y1,y2]
  # fast multiexponential summation method is used
  my $strand=shift;
  my ($x1,$x2,$y1,$y2)=@_;
  my ($sum,$Wx,$Wx_eo,$xsum,$xsum_eo);
	foreach $m (0..$strand->sw->C_last) {
    #initialize Wry(m):
	  my $Wry=$strand->sw->C1($m) * $strand->sw->C2($m) ** ($y2-$x1-2); # m-component of p->omega($y2-$x1)
		#initialize Wx(m):
	  ($Wx,$Wx_eo)=$strand->rev->sum_pceo(
			$strand->revloc($y1,$y2),
			sub {
				my $return=$Wry * $strand->rev->Z10($_[0]+1);
			  $Wry/=$strand->sw->C2($m); # preparing for next ry
				$return;
			},
		);
		# sum over x:
	  ($xsum,$xsum_eo)=$strand->sum_pceo(
      $x1,$x2,
			sub {
				my $return=$strand->Z10($_[0]+1) * $Wx;
			  $Wx/=$strand->sw->C2($m); # preparing for next x
				$return;
			},
		);
		# add to total:
		$sum+=$xsum;
	}		
  $strand->divide_by_Ztotal($sum, $Wx_eo+$xsum_eo);
}
sub upper_limit_sum_p_bubble {
  my $strand=shift;
  my ($x1,$x2,$y1,$y2)=@_;
  my ($sx,$eox)=$strand->sum_loc_Z10_cache($x1,$x2);
  my ($sy,$eoy)=$strand->rev->sum_loc_Z10_cache( $strand->revloc($y1,$y2) );
  $strand->divide_by_Ztotal( $sx * $strand->sw->omega($y1-$x2) * $sy, $eox+$eoy);
}
sub lower_limit_sum_p_bubble {
  my $strand=shift;
  my ($x1,$x2,$y1,$y2)=@_;
  my ($sx,$eox)=$strand->sum_loc_Z10_cache($x1,$x2);
  my ($sy,$eoy)=$strand->rev->sum_loc_Z10_cache( $strand->revloc($y1,$y2) );
  $strand->divide_by_Ztotal( $sx * $strand->sw->omega($y2-$x1) * $sy, $eox+$eoy);
}


1;


