#########################
#
#    Bubblyhelix - DBM_object.pm
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

# DBM files are represented as objects with methods that assume that
# the DBM keys can be sorted numerically (e.g. positions on a chromosome)
# and that the DBM values are strings containing tab-separated lists

package DBM_object;

use File::Path 'mkpath';
use DB_File; # use on Mac


sub open { # this is the constructor
  my $class=shift;
  die unless @_;
  my $dbm_name=pop;
  my $path=join "/", @_;
  mkpath([$path],0,0777) if @_;
  my $SELF;
  dbmopen(%$SELF, "$path/$dbm_name", 0644)
    or die "cannot dbmopen $path/$dbm_name: $!";
  bless $SELF, $class;
};
sub close {
  my $SELF=shift;
  dbmclose(%$SELF);
}
sub read {
  my $SELF=shift;
  my ($key)=@_;
  exists $SELF->{$key} ? split "\t", $SELF->{$key} : ();
}
sub write {
  my $SELF=shift;
  my ($key,@values)=@_;
  $SELF->{$key}=join "\t", @values unless (exists $SELF->{$key});
}

sub brackets {
  my $SELF=shift;
  my ($input)=@_;
  my ($x,$y); # smallest interval such that input is in [x,y[ . Can be undef if input outside range
  while (my ($key,$value) = each %$SELF) {
    $x=$key if (!(defined $x) and $key<=$input);
    $y=$key if (!(defined $y) and $key>$input);
    $x=$key if (defined $x and $x<$key and $key<=$input);
    $y=$key if (defined $y and $y>$key and $key>$input);
  }
	($x,$y);
}
sub dumper {
# Example: $Ztot->dumper( sub {print GP "$_[0] ".&sci_not($_[1],$_[2])."\n"} );
  my $SELF=shift;
  my ($callback)=@_;
  foreach $key (sort {$a<=>$b} keys %$SELF) {
    $callback->($key,$SELF->read($key));
  }
}

sub DESTROY { $_[0]->close }

1;


