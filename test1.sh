#########################
#
#    Bubblyhelix - test1.sh
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


# Fastafiles are kept in a folder, e.g. you could create a folder hg_19/
# containing each chromosome sequence. Such a folder must be "reversified",
# meaning that a subfolder is created with fastafiles for the corresponding
# reverse complement sequences. The script tools/reversify.pl does this for
# each fastafile that has not already been reversified.
# Check the content of the fastafiles/revcompl/ folder
# before and after typing the command below:

/usr/bin/perl tools/reversify.pl fastafiles/

# A fastafile and its reverse complement version are given as arguments to
# the main program probabilitymap/probabilitymap.pl. It also gets its
# options from an editable file (probabilitymap/inputoptions_test1.txt).
# And it needs a name of the sequence that is used for naming a new folder
# ("c4781"). The command below computes a probability map for the c4781
# sequence. Check that two new folders are created, one with results and
# another with a partition function database (pfdbm).

/usr/bin/perl probabilitymap/probabilitymap.pl fastafiles/c4781.fa fastafiles/revcompl/c4781.fa c4781 test/inputoptions_test1.txt 

# compute similarly a probability map for the chrMT (mitochondrion)
# sequence, and check the additions to the output folders:

/usr/bin/perl probabilitymap/probabilitymap.pl fastafiles/chrMT.fa fastafiles/revcompl/chrMT.fa chrMT test/inputoptions_test1.txt 

# The script tools/diffdirs.sh finds the difference between the contents of
# two directories. Check that the newly created probability map results are 
# identical to the precalculated results contained in test/facit_test1/:

tools/diffdirs.sh test/results_test1/ test/facit_test1/

# Now let us inspect the contents of the partition function database: pfdbm_0.075/.
# We must convert the binary DBM files to human-readable text files.
# This can be done with the script tools/pfdbm2txt.pl.
# This script dumps an entire pfdbm-directory to a new "copy" directory,
# in which all DBM files are replaced with corresponding .txt files.
# After doing the following command, you may look into the files.

/usr/bin/perl tools/pfdbm2txt.pl pfdbm_0.075/ test/pfdbm_0.075_dump


