/*
#########################
#
#    Bubblyhelix - pfsweep.cpp
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
*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <cmath>


using namespace std;

typedef int position; // position in sequence (1 is usually first)
typedef int indeks; // index for vectors etc. (-ks not -x)
typedef string filename; // file names like "fastafiles/chrMT.fa"
typedef double pfvalue; // the floating point part of a partion function value
typedef int eovalue; // the exponent offset part of a partion function value

int round_et (double x) {
    return (x >= 0 ? floor(x + 0.5) : ceil(x - 0.5)); // round a double to nearest integer
}


int main( int argc, char *argv[] )
{
	// get parameters from command line:
	//
	if (argc < 13) {
		cerr << "\n Usage: ./pfsweep <fasta_file> <param_file> <poslist_file> <start> <next_carry_position> <Zeo> <partial_sum_Z10> <Z01> <Z11> <Z10> <L_tail> <W>\n\n";
		exit( 1 );
	}

	string commandline; // to put argv into a string

	for (indeks i = 1; i < argc; ++i) { // except argv[0]
		commandline += argv[i];
		commandline += " ";
	}

	filename fastafile, parameterfile, stopsfile;
	position start, next_carry;
	eovalue Zeo;
	pfvalue partial_sum_Z10, Z01, Z11, Z10, L_tail;
	
	istringstream argumentsinput( commandline ); // to extract variables from the string
	
	argumentsinput >> fastafile >> parameterfile >> stopsfile
	     >> start >> next_carry >> Zeo >> partial_sum_Z10 >> Z01 >> Z11 >> Z10 >> L_tail;

	// extract trailing arguments into vector W:
	int W_size = argc - 12;
	vector< pfvalue > W( W_size );
	indeks i = 0;
	
	while (argumentsinput >> W[i]) {
		++i;
	}


	
	
	// Read entire sequence from fasta file and put it in a string:
	//
	ifstream fastainput ( fastafile.c_str() ); // open file
	if (!fastainput) {
	    cerr << "File " << fastafile << " could not be opened" << endl;
		exit( 1 );
	}
		
	string line;
	string sequence = "/"; // the sequence (to begin at index=1)
	
	getline(fastainput, line); // skip header line in fastafile	
	while (fastainput >> line) {
		sequence += line;
	}
	fastainput.close();
	
	

	
	// Read parameters from parameterfile:
	//
	string alphabet;
	int alphabet_size;
	int modulus_z;
	int c_size;
	vector< double > c1, c2, s010, s1;
	vector< double > s11; // vector is used as a quick-and-dirty 2D array: s11[i][j] = s11[i * alphabet_size + j]
	
	ifstream parameterinput ( parameterfile.c_str() ); // open file
	if (!parameterinput) {
	    cerr << "File " << parameterfile << " could not be opened" << endl;
		exit( 1 );
	}
		
	string name;
	
	while (parameterinput >> name) {
		if (name == "modulus_z") {
			parameterinput >> modulus_z;
		}
		else if (name == "alphabet") {
			parameterinput >> alphabet_size;
			for (indeks i=0; i < alphabet_size; ++i) {
				char base;
				parameterinput >> base;
				alphabet += base;
			}
		}
		else if (name == "s010") {
			double x;
			for (indeks i = 0; i < alphabet_size; ++i) {
				parameterinput >> x;
				s010.push_back(x);
			}
		}
		else if (name == "s1") {
			double x;
			for (indeks i = 0; i < alphabet_size; ++i) {
				parameterinput >> x;
				s1.push_back(x);
			}
		}
		else if (name == "s11") {
			double x;
			for (indeks i = 0; i < alphabet_size; ++i) {
				parameterinput >> x;
				s11.push_back(x);
			}
		}
		else if (name == "fixman_freire") {
			parameterinput >> c_size;
			if (c_size != W_size) {
				cerr << "c_size != W_size" << endl;
				exit( 1 );
			}
			for (indeks i=0; i < c_size; ++i) {
				double x1, x2;
				string name;
				parameterinput >> name;
				parameterinput >> x1;
				parameterinput >> x2;
				c1.push_back(x1);
				c2.push_back(x2);
			}
		}
	}
	parameterinput.close();
	


	
	// Read stopsfile and put it in a vector:
	//
	ifstream stopinput ( stopsfile.c_str() ); // open file
	if (!stopinput) {
	    cerr << "File " << stopsfile << " could not be opened" << endl;
		exit( 1 );
	}
		
	vector< position > stops; // the stop positions
	position x;
	
	while (stopinput >> x) {
		stops.push_back( x );
	}
	stopinput.close();
	
	
	
	
	// nested iteration over sequence:
	//
	vector< position >::const_iterator nextstop;
	
	// Search the alphabet index of the base:
	char base = static_cast<char>( toupper( sequence[start] ) );
	indeks a = 0;
	while ( base != alphabet[a] ) {
		a++;
	}
	
	for (nextstop = stops.begin(); nextstop != stops.end(); ++nextstop) { // outer loop over stop positions
		for ( position x = start+1; x <= *nextstop; ++x ) { // inner loop over each position
		    
			// Search the alphabet index of the base:
			char base = static_cast<char>( toupper( sequence[x] ) );
			indeks b = 0;
			while ( base != alphabet[b] ) {
			    b++;
			}
			
			// matrix multiplication:
			Z11 = s11[a * alphabet_size + b ]*(s1[a]*Z01+Z11);
			Z01 = L_tail;
			for (indeks m = 0; m < W_size; ++m) {
				Z01 += W[m];
			}
			for (indeks m = 0; m < W_size; ++m) {
			    W[m] = Z10*c1[m]+W[m]*c2[m];
			}
			Z10 = s1[b]*Z11+s010[b]*Z01;
			partial_sum_Z10 += Z10;
			
			// exponent offset:
			if (x == next_carry) {
				int carry = round_et( log10( partial_sum_Z10 ) );
				Zeo += carry;
				double power = pow(10.0, carry);
				partial_sum_Z10 /= power;
				Z01 /= power;
				Z11 /= power;
				Z10 /= power;
				L_tail /= power;
				for (indeks m = 0; m < W_size; ++m) {
					W[m] /= power;
				}
				next_carry += modulus_z;
			}
			a = b;
			
		} // end inner loop
		
		// write a line of output:
		cout << uppercase << setprecision(16) << *nextstop << " " << next_carry << " " << Zeo << " " << partial_sum_Z10
		     << " " << Z01 << " " << Z11 << " " << Z10 << " " << L_tail;
		for (indeks m = 0; m < W_size; ++m) {
			cout << uppercase << setprecision(16) << " " << W[m];
		}
		cout << endl;
		
	    start = *nextstop;
		
	} // end outer loop
	
	return 0;
} // end of main
