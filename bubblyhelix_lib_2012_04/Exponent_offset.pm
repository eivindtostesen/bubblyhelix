#########################
#
#    Bubblyhelix - Exponent_offset.pm
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

# Math for numbers represented as "floating point with exponent offset".
# A pair ($fp,$eo) is a representation of $fp*10**$eo

my $log10=log(10);

sub log10 {
  my ($fp,$eo)=@_;
  $eo+log($fp)/$log10;
}
sub equal_to {
  my ($m1,$e1,$m2,$e2)=@_;
  (($m1/$m2)*10**($e1-$e2) == 1);
}
sub greater_than {
  my ($m1,$e1,$m2,$e2)=@_;
  (($m1/$m2)*10**($e1-$e2) > 1);
}
sub greater_than_or_equal_to {
  my ($m1,$e1,$m2,$e2)=@_;
  (($m1/$m2)*10**($e1-$e2) >= 1);
}
sub less_than {
  my ($m1,$e1,$m2,$e2)=@_;
  (($m1/$m2)*10**($e1-$e2) < 1);
}
sub less_than_or_equal_to {
  my ($m1,$e1,$m2,$e2)=@_;
  (($m1/$m2)*10**($e1-$e2) <= 1);
}
sub sum_eo {
  # EO scheme for summing (unordered) stream of (fp,eo)-pairs
  my ($iterator)=@_;
	my (%sums,@i,$eo,$tot);
	$sums{$i[1]}+=$i[0] while @i=&$iterator;
	foreach (sort {$a<=>$b} keys %sums) {
	  $tot/=10**($_-$eo) if (defined $eo); # don´t the first time
		$tot+=$sums{$_};
		$eo=$_;
	}
	return ($tot,$eo);
}
sub sci_not { # scientific notation
  my ($fp,$eo)=@_;
  $fp=sprintf("%.16e", $fp);
  if ($fp =~ m<([-+\d.,]+)(e|E)([-+\d.,]+)>) {
    my ($mantissa,$eE,$exponent)=($1,$2,$3); # split $fp
    if (exists $cached_expon_sum{$eo}{$exponent}) {
      return "$mantissa"."$eE".$cached_expon_sum{$eo}{$exponent};
    } else {  
      my $exponentsum=$exponent;
      $exponentsum+=$eo if ($mantissa!=0); # add external exponent 
      $exponentsum=sprintf("%0.0f", $exponentsum); # round to nearest integer
      $exponentsum="+"."$exponentsum" if ($exponentsum>=0);  #add + sign
      if ($exponentsum =~ m<([+-])([\d]+)>) {
        my ($sign,$int)=($1,$2);
        $int=sprintf("%02d", $int); # zero padding
        $exponentsum="$sign"."$int",
      }
      $cached_expon_sum{$eo}{$exponent}=$exponentsum;
      return "$mantissa"."$eE"."$exponentsum";   
    }
  }
  warn "sci_not of ($fp,$eo) is undef";
  undef;
}

1;


