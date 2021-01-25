README
======

Bubblyhelix is software for computing DNA melting profiles.

Bubblyhelix can be downloaded from its homepage: http://bubblyhelix.org

This is a beta version (0.9.4). Later versions are expected to come with less bugs, more documentation and better functionality. Please update accordingly (from http://bubblyhelix.org). Feedback from you as a user may stimulate development.

The Bubblyhelix software was designed to work on large genomic DNA sequences, such as whole human chromosomes, using only a typical laptop or similar. The problem with large sequences is the limited amount of RAM memory. Bubblyhelix uses at most a few gigabytes RAM, without trade-offs in numerical or physical accuracy, it is a lossless algorithm.

Bubblyhelix uses Perl, C++, shell script and DBM files. (And in later versions: gnuplot for graphics.)

Bubblyhelix has a command line interface. It is intended for Unix-like systems. It has only been tested on Mac OS X and CentOS (Linux).


FILES
------

Along with this Readme.txt, you have received these files:

-- License.txt:
Describes what you can and cannot do. Bubblyhelix is Copyright 2012 Eivind Tostesen and licensed under GPL version 3.

-- Citation.txt:
How to cite Bubblyhelix, for example when using results in scientific publications.

-- Authors.txt:
Lists contributors and contact information.

-- ChangeLog.txt:
A formatted list of changes.



Installing?
-----------

After unpacking the Bubblyhelix tarball, you only need to compile the C++ program c_code/pfsweep.cpp. The resulting executable file should be called c_code/pfsweep.

The following is an example of how it could be done. The GNU compiler is used and then pfsweep is tested by executing it with no input arguments, to see if it gives an error message about usage:


$ cd c_code/
$ g++ pfsweep.cpp -o pfsweep
$ ./pfsweep

 Usage: ./pfsweep <fasta_file> <param_file> <poslist_file> <start> <next_carry_position> <Zeo> <partial_sum_Z10> <Z01> <Z11> <Z10> <L_tail> <W>

$ cd ..



