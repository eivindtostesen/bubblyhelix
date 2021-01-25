#########################
#
#    Bubblyhelix - Melting_strand_with_DBM.pm
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

package Melting_strand_with_DBM;

@ISA=qw(Melting_strand_region);
use Melting_strand_region;

sub Ztotal_dbm { # get/set
  my $strand=shift;
  $strand->{"Ztotal DBM"}=shift if (@_);
  $strand->{"Ztotal DBM"};
}

sub pfdbm { $_[0]->{"PF DBM"} }

sub initialize_pfdbm {
  # sets {"PF DBM"} and adds pf_first_xpfvector if it is newly created
  my $strand=shift;
  $strand->{"PF DBM"}=shift;
  $strand->pfdbm->write($strand->pf_first_xpfvector);
}
sub pfdbm_writing_positions { # get/set
  my $strand=shift;
  $strand->{"PF DBM writing positions"}=[sort {$a <=> $b} @_] if (@_);
  @{ $strand->{"PF DBM writing positions"} };
}
sub forward_and_reverse {
  my $strand=shift;
  my ($start,$end)=@_;
  
  # fill arrays both ways in the region:
  my @xpfvector1=$strand->pf_fill_arrays($start,$end);
  my @xpfvector2=$strand->rev->pf_fill_arrays( $strand->revloc($start,$end) ); 

  # compute Ztotal:
  my $TC=$strand->sw->TC;
  if (exists $strand->Ztotal_dbm->{$TC}) {
    # then both strands get Ztotal from forward Ztotal_dbm
    $strand->set_Ztotal($strand->Ztotal_dbm->read($TC));
    $strand->rev->set_Ztotal($strand->Ztotal_dbm->read($TC));
    # and if possible, reverse Ztotal_dbm gets Ztotal from the reverse xpfvector2:
    if ($start==1) {
      $strand->rev->Ztotal_dbm->write($TC,$strand->rev->compute_Ztotal(@xpfvector2));
    }
  } elsif (exists $strand->rev->Ztotal_dbm->{$TC}) {
    # then both strands get Ztotal from reverse Ztotal_dbm
    $strand->set_Ztotal($strand->rev->Ztotal_dbm->read($TC));
    $strand->rev->set_Ztotal($strand->rev->Ztotal_dbm->read($TC));
    # and if possible, forward Ztotal_dbm gets Ztotal from the forward xpfvector1:
    if ($end==$strand->length) {
      $strand->Ztotal_dbm->write($TC,$strand->compute_Ztotal(@xpfvector1));
    }
  } elsif ($start==1 and $end==$strand->length) {
    # then both strands get Ztotal from their own xpfvector
    $strand->set_Ztotal($strand->compute_Ztotal(@xpfvector1));
    $strand->rev->set_Ztotal($strand->rev->compute_Ztotal(@xpfvector2));
    # and copy to both Ztotal_dbm:
    $strand->Ztotal_dbm->write($TC,$strand->Ztotal);
    $strand->rev->Ztotal_dbm->write($TC,$strand->rev->Ztotal);
  } elsif ($start-1 >= $strand->length-$end) {
    # then both strands get Ztotal from the forward xpfvector1
    if ($end<$strand->length) {
      # post sweep from $end to strand length:
	    @xpfvector1=$strand->pf_recursion_no_arrays($strand->length, @xpfvector1);
    }  
    $strand->set_Ztotal($strand->compute_Ztotal(@xpfvector1));
    # copy Ztotal to reverse strand:
    $strand->rev->set_Ztotal($strand->Ztotal);
    # and copy to forward Ztotal_dbm:
    $strand->Ztotal_dbm->write($TC,$strand->Ztotal);
  } else {
    # then both strands get Ztotal from the reverse xpfvector2
    if ($start>1) {
      # post sweep from $end to strand length:
	    @xpfvector2=$strand->rev->pf_recursion_no_arrays($strand->length, @xpfvector2);
    }  
    $strand->rev->set_Ztotal($strand->rev->compute_Ztotal(@xpfvector2));
    # copy Ztotal to forward strand:
    $strand->set_Ztotal($strand->rev->Ztotal);
    # and copy to reverse Ztotal_dbm:
    $strand->rev->Ztotal_dbm->write($TC,$strand->rev->Ztotal);
  }
}
sub pf_fill_arrays {
  my $strand=shift;
  my ($start,$end)=@_;

  # erase previous arrays:
  $strand->tabula_rasa;

  # initialize:
  $strand->{"pf array offset"}=$start;
  my $nearest=($strand->pfdbm->brackets($start))[0];
  die if !defined $nearest;
  $strand->sw->set_iterator($nearest); # rewind
	my @xpfvector=($nearest,$strand->pfdbm->read($nearest));

  if ($start>$nearest) {
    # pre sweep from $nearest to $start:
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
sub pf_recursion {
  my $strand=shift;
  my ($callback1,$callback2,$end,@xpfvector)=@_;
	my $start=$xpfvector[0];
  foreach $stop (grep {$start<$_ and $_<=$end and !exists $strand->pfdbm->{$_}} $strand->pfdbm_writing_positions) {
    @xpfvector=$strand->SUPER::pf_recursion($callback1,$callback2,$stop,@xpfvector);
		$strand->pfdbm->write(@xpfvector);
  }
  @xpfvector=$strand->SUPER::pf_recursion($callback1,$callback2,$end,@xpfvector) if ($xpfvector[0]<$end);
  return @xpfvector;
}

1;


