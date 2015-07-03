// Copyright Steven Riley (sr@stevenriley.net)

// This file is part of the library idsource.

// idsource is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This work is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this idsource.  If not, see <http://www.gnu.org/licenses/>.

#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<iostream>
#include<vector>
#include<cmath>

#include<gsl/gsl_matrix.h>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>

#include"utility.h"

using namespace std;

class Parameters {
private:
	int intNumberOfStringParameters;
	int intNumberOfScalarParameters;
	int intNumberOfVectorParameters;
	int intTotalNumberOfParameters;
	string ** ptptParameterNames;
	int * ptIntSizeOfVectorParameters;
	double ** ptptValuesOfVectorParameters;
	double * ptValuesOfScalarParameters;
	int ** ptptLengthOfVectorParameters;
	int * ptTypeSpecificIndex;
	int * ptTypesOfParameters;
	string ** ptptValuesOfStringParameters;
	int match(string s);
public:
	Parameters() {srerror("No default constructor for Parameter");};
	Parameters(string filename);
	~Parameters();
	string get_string_val(string p);
	string * get_string_pt(string p);
	void set_string(string p, string v);
	double get_scalar_val(string p);
	double * get_scalar_pt(string p);
	void set_scalar(string p, double v);
	double get_vector_val(string p, int i);
	double * get_vector_pt(string p);
	void set_vector(string p, int i, double v);
};
