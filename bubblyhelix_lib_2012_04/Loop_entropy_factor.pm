#########################
#
#    Bubblyhelix - Loop_entropy_factor.pm
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


# See also http://stitchprofiles.uio.no/loopentropy.php


package Loop_entropy_factor;

sub I {
# size minus one of multiexponential
  my ($N)=@_;
  return int(1+log(2*$N));
}

sub sequence_length_range {
  my ($N)=@_;
  my $I=&I($N);
  my $maxN=int(exp($I)/2);
  my $minN=1+int(exp($I-1)/2);
  return ($minN,$maxN);
}
sub arrays {
  my ($N, $d, $sigma, $alfa)=@_;

  #calculate A and B arrays:
  my ($A,$B)=&recursive_approximation($N, $alfa);

  #calculate C1 and C2 arrays:
  my ($C1,$C2)=&AB_to_C1C2($A, $B, $d, $sigma);
  
  return ($A,$B,$C1,$C2);
}

sub recursive_approximation {
  my ($N, $alfa)=@_;
  my (@A,@B);
  my $I=&I($N);
  foreach $n (0..$I) {
    $B[$n]=exp($n-$I);
    $A[$n]=exp(1-$alfa*($I-$n));
    foreach $m (0..$n-1) {
      $A[$n]-=$A[$m]*exp(1-exp($m-$n));
    }
  }
  my $k=0;
  foreach $m (0..$I) {
    $k+=$A[$m]*exp(-$B[$m]*exp(3));
  }
  foreach $n (0..$I) {
    $A[$n]/=$k*exp($alfa*3);
  }
  (\@A,\@B);
}

sub algebraic_approximation {
  my ($N, $alfa)=@_;
  my (@A,@B);
  my $I=&I($N);
  foreach $n (0..$I) {
    $B[$n]=exp($n-$I);
    $A[$n]=exp(1-$alfa*($I-$n));
    foreach $j (0..$n-2) {
      $A[$n]-=(1-exp(1))**($n-$j-1)*exp(2-$alpha*($I-$j));
    }
  }
  my $k=0;
  foreach $m (0..$I) {
    $k+=$A[$m]*exp(-$B[$m]*exp(3));
  }
  foreach $n (0..$I) {
    $A[$n]/=$k*exp($alfa*3);
  }
  (\@A,\@B);
}

sub AB_to_C1C2 {
  my ($A, $B, $d, $sigma)=@_;
  my (@C1,@C2);
  my $I=$#$A;
  foreach $m (0..$I) {
    $C1[$m]=$sigma*$$A[$m]*exp(-$d*$$B[$m])*exp(-4*$$B[$m]);
    $C2[$m]=exp(-2*$$B[$m]);
  }
  (\@C1,\@C2);
}

sub omega_C {
  my ($C1,$C2,$x)=@_;
  my $omega;
  foreach $m (0..$#$C1) {
	  $omega+=$C1->[$m] * $C2->[$m] ** ($x-2);
	}
  $omega;    
}


1;