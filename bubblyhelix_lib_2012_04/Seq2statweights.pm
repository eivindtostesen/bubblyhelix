#########################
#
#    Bubblyhelix - Seq2statweights.pm
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

package Seq2statweights;

use Loop_entropy_factor;
use Iupac;

my $bignumber=1e150;
my $smallnumber=1e-290;
$R=1.987;  #the gas constant

##### initialization

sub new {
  my $class=shift @_;
  bless { @_ }, $class;
};

sub empty_cache { $_[0]->{"cache"}=(); }

sub initialization {
  my $self=shift;
  my ($therm_param,$fastafile)=@_;
  require $therm_param; # a module that contains package: Thermodynamic_parameters
  $self->set_seq( Entire_DNA_sequence->new ); # new seq object
  $self->seq->read_file($fastafile);
  $self->initialize_entropic_factors($self->seq->length);
}

### Sequence

sub seq { $_[0]->{"sequence"} }

sub set_seq { $_[0]->{"sequence"}=$_[1] }

sub set_iterator { $_[0]->seq->set_iterator($_[1]) }

### Nearest neighbor thermodynamics

sub TC { $_[0]->{"temperature_celcius"} }

sub TK { $_[0]->{"temperature_kelvin"} }

sub salt1 { $_[0]->{"salt_molar"} }

sub modulus_Z { $_[0]->{"modulus_z"} }

sub s1 { $_[0]->{"%s1"}->{ $_[0]->seq->base($_[1]) } }

sub s010 { $_[0]->{"%s010"}->{ $_[0]->seq->base($_[1]) } }

sub s11 {
  my ($a,$b)=$_[0]->seq->NN($_[1]);
  $_[0]->{"%s11"}->{$a}->{$b};
}

sub next_transfer_matrix {
  my ($a,$b)=$_[0]->{"sequence"}->next_NN;
  return (
    $_[0]->{"%s1"}->{$a},
    $_[0]->{"%s1"}->{$b},
    $_[0]->{"%s11"}->{$a}->{$b},
    $_[0]->{"%s010"}->{$b},
  );
}

sub set_stabilities {
  my $self=shift;
  ($self->{"temperature_celcius"},$self->{"salt_molar"})=@_;
  $self->{"temperature_kelvin"}=273.15+$self->{"temperature_celcius"};
  ($self->{"%s11"},$self->{"%s1"},$self->{"%s010"})=&Thermodynamic_parameters::nested_hashes($self->{"temperature_celcius"},$self->{"salt_molar"});
  &Iupac::include($self->{"%s11"},$self->{"%s1"},$self->{"%s010"},[ $self->seq->incompletely_specified ]);
  my $s_GCGC=sqrt($self->{"%s11"}->{"G"}->{"C"}*$self->{"%s11"}->{"C"}->{"G"});
  #my $s_TATA=sqrt($self->{"%p11"}->{"T"}->{"A"}*$self->{"%p11"}->{"A"}->{"T"});
  $self->{"modulus_z"}=sprintf("%d", 1+log($bignumber)/log($s_GCGC));
  die "Extreme temperature!" if ($self->{"modulus_z"}<5);
}


### Loop entropy factor (and tail entropy)

sub initialize_entropic_factors {
  my $self=shift;
  my ($N)=@_; # $N is length
  $self->empty_cache;
  ($self->{"d"},$self->{"sigma"},$self->{"alpha"})=&Thermodynamic_parameters::loop_parameters();
  ($self->{"sigma_tail"})=&Thermodynamic_parameters::tail_parameters();
  ($self->{"A"},$self->{"B"},$self->{"C1"},$self->{"C2"})=&Loop_entropy_factor::arrays($N,$self->{"d"},$self->{"sigma"},$self->{"alpha"});
  $self->{"sequence_length_range"}=[ &Loop_entropy_factor::sequence_length_range($N) ];
}


sub sigma_tail { $_[0]->{"sigma_tail"} }

sub C_last { $#{$_[0]->{"C1"}} } # index of last element

sub C1 { $_[0]->{"C1"}->[$_[1]] }

sub C2 { $_[0]->{"C2"}->[$_[1]] }

sub A { $_[0]->{"A"}->[$_[1]] }

sub B { $_[0]->{"B"}->[$_[1]] }

sub omega {
  $_[0]->{"cache"}->{"omega"}->{$_[1]}=&Loop_entropy_factor::omega_C($_[0]->{"C1"},$_[0]->{"C2"},$_[1]) unless (exists $_[0]->{"cache"}->{"omega"}->{$_[1]});
  $_[0]->{"cache"}->{"omega"}->{$_[1]};    
}


### IO

sub parameter_file { $_[0]->{"parameter file"} }


sub read_parameter_file_v3 {
  my $self=shift;
  ($self->{"parameter file"})=@_;

  my (%base_index,@s11,@s1,@s010,$alphabet_size,@alphabet);
  
  open(P,$self->{"parameter file"}) ||
    die("Cannot open file \'$_[0]\' ");
  
  # erase previous:
  delete $self->{"C1"};
  delete $self->{"C2"};
  
  while(my $line=<P>){
    chomp($line);
    next if $line eq "";
    my ($name,@values)=split " ", $line;
    ($self->{"temperature_celcius"})=@values if ($name eq "temperature_celcius");
    ($self->{"salt_molar"})=@values if ($name eq "salt_molar");
    ($self->{"modulus_z"})=@values if ($name eq "modulus_z");
    ($self->{"alpha"})=@values if ($name eq "alpha");
    ($self->{"sigma"})=@values if ($name eq "sigma");
    ($self->{"d"})=@values if ($name eq "d");
    ($self->{"sigma_tail"})=@values if ($name eq "sigma_tail" or $name eq "s_tail");
    push (@s11, \@values) if ($name eq "s11");
    @s1=@values if ($name eq "s1");
    @s010=@values if ($name eq "s010");
    ($self->{"sequence_length_range"})=@values if ($name eq "sequence_length_range");
    if($name eq "c1c2"){
      push (@{ $self->{"C1"} }, $values[0]);
      push (@{ $self->{"C2"} }, $values[1]);
    }
    if($name eq "alphabet"){
      ($alphabet_size,@alphabet)=@values;
      %base_index=map {$alphabet[$_],$_} 0..$#alphabet; # i.e. (A,0,C,1,G,2,T,3,N,4)
    }
  }
  close P;
  
  # from arrays to nested hashes:
  $self->{"%s1"}={ map {$alphabet[$_],$s1[$_]} 0..$#alphabet };
  $self->{"%s010"}={ map {$alphabet[$_],$s010[$_]} 0..$#alphabet };
  foreach $i (0..$#alphabet) {
    foreach $j (0..$#alphabet) {
      $self->{"%s11"}->{$alphabet[$i]}->{$alphabet[$j]}=$s11[$i][$j];
    }
  }
}


sub write_parameter_file_v3 {  
  my $self=shift;
  ($self->{"parameter file"})=@_;

  open my $p,">",$self->{"parameter file"}
    or die "cannot open file: $!";

  print {$p} "poland_scheraga parameters_v3\n\n";
  ############    NN Stabilities:
  print {$p} "temperature_celcius ".$self->TC."\n";
  print {$p} "salt_molar ".$self->salt1."\n";
  print {$p} "modulus_z ".$self->modulus_Z."\n";
  my @alphabet=qw/A C G T/;
  push @alphabet, grep {/[^acgtACGT]/} keys %{$self->{"%s11"}};
  print {$p} "alphabet ".@alphabet." ".(join " ", @alphabet)."\n\n";
  foreach $X (@alphabet) {
    print {$p} "s11";
    foreach $Y (@alphabet) {
      print {$p} " ".$self->{"%s11"}->{$X}->{$Y};
    }
    print {$p} "\n";
  }
  print {$p} "\n";
  print {$p} "s010 ".(join " ", map {$self->{"%s010"}->{$_}} @alphabet );
  print {$p} "\n";
  print {$p} "\n";
  print {$p} "s1 ".(join " ", map {$self->{"%s1"}->{$_}} @alphabet );
  print {$p} "\n";
  ############    Loop/tail entropies:
  print {$p} "\nalpha ".$self->{"alpha"}."\n";
  print {$p} "sigma ".$self->{"sigma"}."\n";
  print {$p} "d ".$self->{"d"}."\n";
  print {$p} "s_tail ".$self->{"sigma_tail"}."\n";
  print {$p} "\nsequence_length_range ".(join " ",@{$self->{"sequence_length_range"}})."\n";
  print {$p} "fixman_freire ".@{$self->{"C1"}}."\n";
  foreach $m (0..$#{$self->{"C1"}}) {
    print {$p} "c1c2 ".$self->{"C1"}->[$m]." ".$self->{"C2"}->[$m]."\n";
  }
  close $p;
}

1;

