#########################
#
#    Bubblyhelix - Melting_strand_with_DBM_C.pm
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

package Melting_strand_with_DBM_C;

@ISA=qw(Melting_strand_with_DBM);
use Melting_strand_with_DBM;

# where to find the sweeper program:
my $sweeper="c_code/pfsweep";

sub pf_recursion_no_arrays {
  my $strand=shift;
  my ($end,@xpfvector)=@_;
  # check existence of sweeper program:
  die unless (-e $sweeper);

	my $start=$xpfvector[0];
  my @stops=grep {$start<$_ and $_<=$end and !exists $strand->pfdbm->{$_}} $strand->pfdbm_writing_positions;

  # create temporary file with stop positions to be read by sweeper program:
	my $stops_file="stop_positions_$start-$end.r".int(rand(1000)).".$$.tmp";
	open TMP, ">$stops_file" or
  	die("cannot create $stops_file");
	map {print TMP "$_\n"} @stops;
	print TMP "$end\n" if (!@stops or $stops[-1]!=$end);
	close TMP;

  # the sweeper program needs two more files:
  my $fasta_file=$strand->sw->seq->fastafile;
  my $parameter_file=$strand->sw->parameter_file;

  # run the sweeper program, read its output lines through pipe, and write them to pfdbm:
	open C, "./$sweeper $fasta_file $parameter_file $stops_file @xpfvector |"
		or die "cannot pipe from $sweeper: $!";
	while (my $line=<C>) {
		chomp($line);
		@xpfvector=split " ", $line;
		$strand->pfdbm->write(@xpfvector) if (@stops and $xpfvector[0]<=$stops[-1]);
	}
	close C;
  
	unlink $stops_file;
  $strand->sw->set_iterator($end); # bring it forward to position $end
  return @xpfvector;
}

1;


