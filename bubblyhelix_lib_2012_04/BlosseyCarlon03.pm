#########################
#
#    Bubblyhelix - BlosseyCarlon03.pm
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

#This is the parameter set from R. Blossey and E. Carlon (2003): Reparametrizing loop entropy weights: Effect on DNA melting curves. Phys. Rev. E, 68, 061911.
#It's a reparametrization of the parameters of R. D. Blake and S. G. Delcourt (1998): Thermal stability of DNA. Nucl. Acids Res. 26, 3323-3332.
#such that alfa=2.15 and sigma=1.26e-4

#The usage is: open-to-closed: sij=exp((-dHij/RTij)*(1-Tij/T)),
#where Tij=T0ij+log10[Na+]*dTij/dlog10[Na+]. T0ij is for 1.0 M salt.
#Below are tables of dHij, T0ij and dTij/dlog10[Na+] (in kcal/mol and degrees C).


package Thermodynamic_parameters;

%T0ij=(
  "AA=TT" => 89.08,
  "AT" => 89.38,
  "TA" => 79.47,
  "CA=TG" => 89.71,
  "AC=GT" => 121.17,
  "AG=CT" => 98.49,
  "GA=TC" => 105.09,
  "CG" => 105.28,
  "GC" => 143.73,
  "CC=GG" => 118.49,
);

%dHij=(
  "AC=GT" => 10.51,
  "TA" => 7.81,
  "AA=TT" => 8.45,
  "GC" => 11.91,
  "CC=GG" => 10.34,
  "GA=TC" => 9.47,
  "AT" => 8.50,
  "CG" => 9.53,
  "AG=CT" => 9.10,
  "CA=TG" => 8.51,
);

%dTij_dlog10Na=(
  "AC=GT" => 13.71,
  "TA" => 22.08,
  "AA=TT" => 19.78,
  "GC" => 10.21,
  "CC=GG" => 14.18,
  "GA=TC" => 16.94,
  "AT" => 19.19,
  "CG" => 16.01,
  "AG=CT" => 17.30,
  "CA=TG" => 19.35,
);

%nn2ten=(
  AA => "AA=TT",
  TT => "AA=TT",
  AT => "AT",
  TA => "TA",
  CA => "CA=TG",
  TG => "CA=TG",
  "GT" => "AC=GT",
  AC => "AC=GT",
  CT => "AG=CT",
  AG => "AG=CT",
  GA => "GA=TC",
  TC => "GA=TC",
  CG => "CG",
  GC => "GC",
  GG => "CC=GG",
  CC => "CC=GG",
);

$R=1.987;  #the gas constant

my @four=qw/A C G T/;
my @sixteen=qw/AA TT AT TA CA TG GT AC CT AG GA TC CG GC GG CC/;

($d, $sigma, $alfa)=(1,1.26e-4,2.15); # loop entropy parameters

$s_tail="1.0"; # tail entropy parameter

sub loop_parameters {
  ($d, $sigma, $alfa);
}  
sub tail_parameters {
  ($s_tail);
}  

sub statistical_weight_hashes {
# Calculation of salt- and temperature dependent statistical weights
  my ($TC,$Na_conc)=@_;
  my (%s_11,%s_1,%s_010);
 
  %s_1=map {$_,1} @four;
  %s_010=map {$_,0} @four;
  
  while (($key,$value)=each %T0ij) {
    $Tij{$key}=$value+273.15+(log($Na_conc)/log(10))*$dTij_dlog10Na{$key};
  }
  my $T=$TC+273.15;
  foreach $NN (@sixteen) {
    my $ij=$nn2ten{$NN};
    $s_11{$NN}=exp((-1000*$dHij{$ij}/($R*$Tij{$ij}))*(1-$Tij{$ij}/$T));
  }
  return (\%s_11,\%s_1,\%s_010);
}

sub nested_hashes {
# Calculation of salt- and temperature dependent statistical weights
  my ($TC,$Na_conc)=@_;
  my (%s_11,%s_1,%s_010);
  
  %s_1=map {$_,1} @four;
  %s_010=map {$_,0} @four;
  
  while (($key,$value)=each %T0ij) {
    $Tij{$key}=$value+273.15+(log($Na_conc)/log(10))*$dTij_dlog10Na{$key};
  }
  my $T=$TC+273.15;
  foreach $X (@four) {
    foreach $Y (@four) {
      my $ij=$nn2ten{$X.$Y};
      $s_11{$X}->{$Y}=exp((-1000*$dHij{$ij}/($R*$Tij{$ij}))*(1-$Tij{$ij}/$T));
    }  
  }
  return (\%s_11,\%s_1,\%s_010);
}


sub nested_dHdS_hashes {
# Calculation of salt-dependent dH11 and dS11
  my ($Na_conc)=@_;
  my (%dH_11,%dS_11);
  
  while (($key,$value)=each %T0ij) {
    $Tij{$key}=$value+273.15+(log($Na_conc)/log(10))*$dTij_dlog10Na{$key};
  }
  foreach $X (@four) {
    foreach $Y (@four) {
      my $ij=$nn2ten{$X.$Y};
      $dH_11{$X}->{$Y}=-1000*$dHij{$ij};
      $dS_11{$X}->{$Y}=-1000*$dHij{$ij}/$Tij{$ij};
    }  
  }
  return (\%dH_11,\%dS_11);
}


1;

