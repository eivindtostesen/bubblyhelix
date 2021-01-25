#########################
#
#    Bubblyhelix - Entire_DNA_sequence.pm
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

package Entire_DNA_sequence;

@ISA=qw(Antiparallel_pair);
use Antiparallel_pair;


sub new {
  my $class=shift @_;
  bless { @_ }, $class;
};
sub read_file {
  my $seq=shift;
  $seq->fastafile(shift) if (@_);
  open(SEQFIL, $seq->{"fastafile"}) ||
    die("Cannot open file");
  my $line=<SEQFIL>;  #use if first line is a header
  while($line=<SEQFIL>){
    chomp($line);
	$seq->{"string"}.=uc($line);
  }
  close SEQFIL;
  $seq->length(length($seq->{"string"}));
}
sub fastafile { # get/set
  my $seq=shift;
  $seq->{"fastafile"}=shift if (@_);
  $seq->{"fastafile"};
}
sub length { # get/set
  my $seq=shift;
  $seq->{"length"}=shift if (@_);
  $seq->{"length"};
}
sub whole_sequence {
  my $seq=shift;
  ($seq->{"string"},$seq->{"length"});
}
sub sub_seq {
  my $seq=shift;
  my ($start,$end)=@_;
  substr($seq->{"string"},$start-1,$end-$start+1);
}

sub base { substr($_[0]->{"string"},$_[1]-1,1) }

sub NN { split "", substr($_[0]->{"string"},$_[1]-2,2) }

sub next_base { $_[0]->base(++$_[0]->{"_iter"}) } # preincrement the iterator

sub next_NN { $_[0]->NN(++$_[0]->{"_iter"}) } # preincrement the iterator

sub set_iterator {
  my $seq=shift;
  my ($x)=@_;
  $seq->{"_iter"}=$x;
}
sub first_base {
  my $seq=shift;
  $seq->base(1);
}
sub last_base {
  my $seq=shift;
  $seq->base($seq->length);
}

########  Sequence statistics

sub n_flanks {
  my $seq=shift;
  if (!exists $seq->{"stats"}->{"n_flanks"}) {
    my ($left,$right)=(0,0);
    $left++ while ("N" eq $seq->base($left+1));
    $right++ while ("N" eq $seq->base($seq->length-$right));
    $seq->{"stats"}->{"n_flanks"}=[ $left,$right ];
  }
  return @{ $seq->{"stats"}->{"n_flanks"} };
}
sub stats { # collects statistics (optionally in region $start..$end)
  my $seq=shift;
  my ($start,$end)=@_ if (@_);
  $start=1 if (!$start or $start<1);
  $end=$seq->length if (!$end or $end>$seq->length);
  if (!exists $seq->{"stats"}->{"$start $end"}) {
    my ($a,$b);
    $a=$seq->base($start);
    $seq->{"stats"}->{"$start $end"}->{"first base"}=$a;
    $seq->{"stats"}->{"$start $end"}->{"symbol counts"}->{$a}++;
    foreach ($start+1..$end) {
      $b=substr($seq->{"string"},$_-1,1);
      $seq->{"stats"}->{"$start $end"}->{"symbol counts"}->{$b}++;
      $seq->{"stats"}->{"$start $end"}->{"NN counts"}->{$a}->{$b}++;
      $a=$b;
    }
    $seq->{"stats"}->{"$start $end"}->{"last base"}=$b;
  }
  return $seq->{"stats"}->{"$start $end"};
}
sub symbol_counts { # base composition (optionally in region $start..$end)
  my $seq=shift;
  my $s=$seq->stats(@_)->{"symbol counts"};
  return map {$_,$s->{$_}} sort {$s->{$b} <=> $s->{$a}} keys %$s;
}
sub nn_counts { # doublet composition (optionally in region $start..$end)
  my $seq=shift;
  my $nn=$seq->stats(@_)->{"NN counts"};
  my %nnc= map { $_a=$_; map {$_a.$_, $nn->{$_a}->{$_} } keys %{ $nn->{$_a} } } keys %$nn;
  return map {$_,$nnc{$_}} sort {$nnc{$b} <=> $nnc{$a}} keys %nnc;
  #return map { $a=$_; map {$a.$_, $nn->{$a}->{$_} } sort keys %{ $nn->{$a} } } sort keys %$nn;
}
sub symbols { # alphabet (optionally in region $start..$end)
  my $seq=shift;
  sort keys %{$seq->stats(@_)->{"symbol counts"}};
}
sub incompletely_specified { # uncertain bases (optionally in region $start..$end)
  my $seq=shift;
  sort grep {/[^acgtACGT]/} $seq->symbols(@_);
}
sub gc_content { # (optionally in region $start..$end)
  my $seq=shift;
  my $s=$seq->stats(@_)->{"symbol counts"};
  my ($total,$gc);
  map {$gc+=$s->{$_}} grep {/[cgCG]/} keys %{$s};
  map {$total+=$s->{$_}} grep {/[acgtACGT]/} keys %{$s};
  return $gc/$total unless ($total==0);
}


1;


