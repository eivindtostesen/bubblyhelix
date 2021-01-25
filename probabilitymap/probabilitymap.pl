#!/usr/bin/perl -w

#########################
#
#    Bubblyhelix - probabilitymap.pl
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

# usage: perl probabilitymap/probabilitymap.pl <forw.fa> <rev.fa> <seqname> <optionsfile>


use FindBin qw($Bin);
use lib "$Bin/../bubblyhelix_lib_2012_04";

#use Melting_strand;
#use Melting_strand_region;
#use Melting_strand_with_DBM;
use Melting_strand_with_DBM_C;

use Entire_DNA_sequence;
use Seq2statweights;
use Exponent_offset;
use DBM_object;
use File::Path 'mkpath';

my %module_file=(
  blossey_and_carlon_2003 => "BlosseyCarlon03.pm",
  blake_and_delcourt_1998 => "BlakeDelcourt98.pm",
);


############################################ OPTIONS/ARGUMENTS

my ($forward_fa,$reverse_fa,$chr,$options_file)=@ARGV;

# READ OPTIONS FILE:

my ($pfdb,$resdb,$pfdbm_divisor,$lines_per_file,$pformat,@temps);

open(OPT,$options_file) ||
  die("Cannot open file \'$options_file\' ");
while(my $line=<OPT>){
  chomp($line);
  next if $line eq "";
  my ($name,@values)=split " ", $line;
  ($therm_param)=@values if ($name eq "thermodynamic_parameters");
  ($Na_conc)=@values if ($name eq "salt_molar");
  ($pfdb)=@values if ($name eq "dirname_pfdbm");
  ($resdb)=@values if ($name eq "dirname_results");
  ($pfdbm_divisor)=@values if ($name eq "pfdbm_divisor");
  ($lines_per_file)=@values if ($name eq "output_lines_per_file");
  ($pformat)=@values if ($name eq "output_prob_format");
  ($output_postprocess)=@values if ($name eq "output_postprocess");
  @temperatures=@values if ($name eq "list_temperatures_celcius");
}
close OPT;


#############################################    GET SEQUENCE AND THERMODYNAMIC PARAMETERS

# create pair of melting strands (forward and reverse):
my ($LR,$RL)=Melting_strand_with_DBM_C->new_pair;

# sequence and parameters for $LR:
$LR->set_sw( Seq2statweights->new );
$LR->sw->initialization($module_file{$therm_param},$forward_fa);

# sequence and parameters for $RL:
$RL->set_sw( Seq2statweights->new );
$RL->sw->initialization($module_file{$therm_param},$reverse_fa);

# set length of the pair ($LR,$RL) (after having read the fastafiles):
$LR->set_pair_length($LR->sw->seq->length);

##########

# set pfdbm writing positions:
my $db_interval=sprintf "%.0e", $LR->length/$pfdbm_divisor;
my ($n_left,$n_right)=$LR->sw->seq->n_flanks;
$LR->pfdbm_writing_positions( map {$_*$db_interval} int(1+$n_left/$db_interval)..int(($LR->length-$n_right)/$db_interval) );
$RL->pfdbm_writing_positions( $RL->revloc( $LR->pfdbm_writing_positions ) );

# create folders to put log files, parameters files, DBMs, etc.:
mkpath(
  [join "/",$Bin,$pfdb,$chr,"parameters"],
  0,
  0777
);


# write to log file:
open LOG, ">>$Bin/$pfdb/$chr/pfdbm_log.txt"
  or die "cannot open file: $!";
print LOG "date ".`date`;
print LOG "fastafile_forward $forward_fa\n";
print LOG "fastafile_reverse $reverse_fa\n";
print LOG "thermodynamic_parameters $therm_param\n";
print LOG "salt_molar $Na_conc\n";
print LOG "sequence_length ".$LR->length."\n";
print LOG "sequence_composition_forward ".(join " ", $LR->sw->seq->symbol_counts)."\n";
print LOG "sequence_composition_reverse ".(join " ", $RL->sw->seq->symbol_counts)."\n";
print LOG "nn_composition_forward ".(join " ", $LR->sw->seq->nn_counts)."\n";
print LOG "nn_composition_reverse ".(join " ", $RL->sw->seq->nn_counts)."\n";
print LOG "pfdbm_divisor $pfdbm_divisor\n";
print LOG "writing_positions_forward ".(join " ",$LR->pfdbm_writing_positions)."\n";
print LOG "writing_positions_reverse ".(join " ",$RL->pfdbm_writing_positions)."\n";
print LOG "\n";
close LOG;


# open or create Ztotal DBMs:
$LR->Ztotal_dbm( DBM_object->open($Bin,$pfdb,$chr,"Ztotal_forward") );
$RL->Ztotal_dbm( DBM_object->open($Bin,$pfdb,$chr,"Ztotal_reverse") );

# compute probability profile for each temperature:
foreach $TC ((@temperatures)) {
  
  # read or write parameter file (input to C program):
  my $parfile=join "/",$Bin,$pfdb,$chr,"parameters","$TC.txt";
  if (-f $parfile) {
    # read the existing file:
    $LR->sw->read_parameter_file_v3($parfile);
    $RL->sw->read_parameter_file_v3($parfile);
  } else {
    # update temperature:
    $LR->sw->set_stabilities($TC,$Na_conc);
    $RL->sw->set_stabilities($TC,$Na_conc);
    # create a new file:
    $LR->sw->write_parameter_file_v3($parfile);
    $RL->sw->write_parameter_file_v3($parfile);
  }
  
  # open or create pfdbm:
  $LR->initialize_pfdbm( DBM_object->open($Bin,$pfdb,$chr,$TC,"forward") );
  $RL->initialize_pfdbm( DBM_object->open($Bin,$pfdb,$chr,$TC,"reverse") );

  # create folder to put the output files:
  mkpath(
    [(join "/",$Bin,$resdb,"probabilities",$TC)],
    0,
    0777
  );
  
  # create output file for each sequence window [start,end]:
  foreach $filenumb (reverse(1..int($LR->length/$lines_per_file)+1)) {
    my ($start,$end)=(($filenumb-1)*$lines_per_file,$filenumb*$lines_per_file-1);
    $start=1 if ($start==0);
    $end=$LR->length if ($end>$LR->length);
    my $outputfile=$chr.".".$start."_".$end.".txt";
    my $path=join "/",$Bin,$resdb,"probabilities",$TC,$outputfile;
    next if (-e $path or -e "$path.gz");
    open PP, ">$path"
      or die "cannot create file: $!";
    $LR->forward_and_reverse($start,$end); # compute partition functions
    $LR->p_1_profile($start,$end,sub {printf PP "$pformat\n", &sci_not(@_)});
    close PP;
	system "$output_postprocess $path" if (defined $output_postprocess and $output_postprocess eq "gzip");
  }
  
  $LR->pfdbm->close;
  $RL->pfdbm->close;
}

## the end
