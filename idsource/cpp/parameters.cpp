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

#include"parameters.h"
#include"utility.h"

// A generic parameter class to read froma  file and give pointers.
// Includes vector parameters
Parameters::Parameters(string filename) {

	// Read in lines from filename and process them as parameters
	// Parameter types are: 0 - string, 1 - scalar, 2 - vector
	// This is the only constructor for this class
	// The memory structure of an instance of the class cannot be changed after it has been loaded

	// Initialize housekeeping variables
	intNumberOfStringParameters = 0;
	intNumberOfScalarParameters = 0;
	intNumberOfVectorParameters = 0;
	intTotalNumberOfParameters = 0;

	// Declare some local variables
	string line,name,stringval,junk,strtmp,strtmp2,filecontents;
	int ptype;
	char ch;

	// Open the file and put the contents into a single string
	// Ignore anything within square brackets
	ifstream ifs(filename.c_str());
	if (ifs.fail()) srerror("SREC:001 Can't open the file in constructor for PopulationList: "+filename);
	filecontents = "";
	while (ifs) {
		ifs.get(ch);
		if (ch != '[') filecontents = filecontents + ch;
		else {
			ifs.get(ch);
			while (ch != ']') ifs.get(ch);
		}
	}
	ifs.close();

	// Read each line of the file contents and parse it
	stringstream ifs1(filecontents);
	while (getline(ifs1,line)) {

		// Break up the line
		stringstream ss(line);
		if (line != "") {

			getline(ss,name,',');
			getline(ss,strtmp,',');
			stringstream sstype(strtmp);
			sstype >> ptype;

			if (ptype==0) intNumberOfStringParameters++;
			else if (ptype==1) intNumberOfScalarParameters++;
			else if (ptype==2) intNumberOfVectorParameters++;
			else srerror("Parameter type must be either 0, 1, or 2");

		}
	}

	// Process the first sweep and initialize the arrays and tmp variables
	intTotalNumberOfParameters = 	intNumberOfStringParameters +
									intNumberOfScalarParameters +
									intNumberOfVectorParameters;

	// Allocate all memory here
	ptptParameterNames = new string * [intTotalNumberOfParameters];
	for (int i=0;i<intTotalNumberOfParameters;i++) ptptParameterNames[i] = new string;
	ptptValuesOfStringParameters = new string * [intNumberOfStringParameters];
	for (int i=0;i<intNumberOfStringParameters;i++) ptptValuesOfStringParameters[i] = new string;
	ptIntSizeOfVectorParameters = new int[intNumberOfVectorParameters];
	ptptValuesOfVectorParameters = new double * [intNumberOfVectorParameters];
	ptValuesOfScalarParameters = new double [intNumberOfScalarParameters];
	ptTypeSpecificIndex = new int [intTotalNumberOfParameters];
	ptTypesOfParameters = new int [intTotalNumberOfParameters];

	// Open the file to read again and process each of the variable types
	int intTmpCountStrings = 0;
	int intTmpCountScalars = 0;
	int intTmpCountVectors = 0;

	// Read each line again and process the values
	stringstream ifs2(filecontents);
	for (int i=0;i<intTotalNumberOfParameters;i++) {

		// Let the number of
		getline(ifs2,line);

		// Break up the line and read the name, type and (first value)
		stringstream ss(line);
		getline(ss,name,',');
		*(ptptParameterNames[i]) = name;
		getline(ss,strtmp,',');
		stringstream sstype(strtmp);
		sstype >> ptype;
		ptTypesOfParameters[i] = ptype;

		if (ptype==0) {
			getline(ss,strtmp,',');
			*(ptptValuesOfStringParameters[intTmpCountStrings]) = strtmp;
			ptTypeSpecificIndex[i] = intTmpCountStrings;
			intTmpCountStrings++;
		}
		else if (ptype==1) {
			// Up to here XX next thing is to be able to handle the numeric parameters
			// next thing is a temp string stream to format the data
			getline(ss,strtmp,',');
			stringstream ssval(strtmp);
			ssval >> ptValuesOfScalarParameters[intTmpCountScalars];
			ptTypeSpecificIndex[i] = intTmpCountScalars;
			intTmpCountScalars++;
		}
		else if (ptype==2) {

			// Count the size of the scalar and then read in its values
			getline(ss,strtmp);
			stringstream ssval1(strtmp),ssval2(strtmp);
			int intSizeCurrent = 0;
			while (getline(ssval1,junk,',')) {
				if (junk != "") intSizeCurrent++;
			}
			ptIntSizeOfVectorParameters[intTmpCountVectors] = intSizeCurrent;
			ptptValuesOfVectorParameters[intTmpCountVectors] = new double [intSizeCurrent];
			int intWhileCounter = 0;
			while (getline(ssval2,strtmp2,',')) {
				if (strtmp2 != "") {
					stringstream ssval3(strtmp2);
					ssval3 >> *(ptptValuesOfVectorParameters[intTmpCountVectors]+intWhileCounter);
					intWhileCounter++;
				}
			}
			ptTypeSpecificIndex[i] = intTmpCountVectors;
			intTmpCountVectors++;
		}
		else srerror("Parameter type must be either 0, 1, or 2");
	}

};

Parameters::~Parameters() {
	for (int i=0;i<intTotalNumberOfParameters;i++) delete ptptParameterNames[i];
	for (int i=0;i<intNumberOfStringParameters;i++) delete ptptValuesOfStringParameters[i];
	for (int i=0;i<intNumberOfVectorParameters;i++) delete [] ptptValuesOfVectorParameters[i];
	delete [] ptptParameterNames;
	delete [] ptptValuesOfStringParameters;
	delete [] ptIntSizeOfVectorParameters;
	delete [] ptptValuesOfVectorParameters;
	delete [] ptValuesOfScalarParameters;
	delete [] ptTypeSpecificIndex;
	delete [] ptTypesOfParameters;
};

int Parameters::match(string s) {
	int rtn = 0;
	bool flagdontstop = 1;
	while (rtn < intTotalNumberOfParameters && flagdontstop) {
		if (*(ptptParameterNames[rtn]) == s) flagdontstop = 0;
		else rtn++;
	}
	if (rtn == intTotalNumberOfParameters) srerror("Attempting to match an unknown parameter string: "+s);
	return rtn;
}

string Parameters::get_string_val(string p) {
	string s;
	int genindex = match(p);
	if (ptTypesOfParameters[genindex]!=0) srerror("get_string_val is only for strings.");
	int stringindex = ptTypeSpecificIndex[genindex];
	s = (*ptptValuesOfStringParameters[stringindex]);
	return s;
};

string * Parameters::get_string_pt(string p) {
	int genindex, stringindex;
	genindex = match(p);
	if (ptTypesOfParameters[genindex]!=0) srerror("get_string_pt is only for strings.");
	stringindex = ptTypeSpecificIndex[genindex];
	return ptptValuesOfStringParameters[stringindex];
};

void Parameters::set_string(string p, string v){
	int genindex, stringindex;
	genindex = match(p);
	if (ptTypesOfParameters[genindex]!=0) srerror("set_string is only for strings.");
	stringindex = ptTypeSpecificIndex[genindex];
	*(ptptValuesOfStringParameters[stringindex]) = v;
};

double Parameters::get_scalar_val(string p) {
	double val;
	int genindex = match(p);
	if (ptTypesOfParameters[genindex]!=1) srerror("get_scalar_val is only for scalars.");
	int scalarindex = ptTypeSpecificIndex[genindex];
	val = ptValuesOfScalarParameters[scalarindex];
	return val;
};

double * Parameters::get_scalar_pt(string p) {
	int genindex, scalarindex;
	genindex = match(p);
	if (ptTypesOfParameters[genindex]!=1) srerror("get_scalar_pt is only for scalars.");
	scalarindex = ptTypeSpecificIndex[genindex];
	return ptValuesOfScalarParameters+scalarindex;
};

void Parameters::set_scalar(string p, double v){
	int genindex, scalarindex;
	genindex = match(p);
	if (ptTypesOfParameters[genindex]!=1) srerror("set_scalar is only for scalars.");
	scalarindex = ptTypeSpecificIndex[genindex];
	ptValuesOfScalarParameters[scalarindex] = v;
};

double Parameters::get_vector_val(string p, int i) {
	double val;
	int genindex = match(p);
	if (ptTypesOfParameters[genindex]!=2) srerror("get_vector_val is only for vectors.");
	int vectorindex = ptTypeSpecificIndex[genindex];
	int veclength = ptIntSizeOfVectorParameters[vectorindex];
	if (i < 0 || i >= veclength) srerror("out of range in get_vector_val.");
	val = (ptptValuesOfVectorParameters[vectorindex])[i];
	return val;
};

double * Parameters::get_vector_pt(string p) {
	double * val;
	int genindex = match(p);
	if (ptTypesOfParameters[genindex]!=2) srerror("get_vector_pt is only for vectors.");
	int vectorindex = ptTypeSpecificIndex[genindex];
	val = ptptValuesOfVectorParameters[vectorindex];
	return val;
};

void Parameters::set_vector(string p, int i, double v) {
	int genindex = match(p);
	if (ptTypesOfParameters[genindex]!=2) srerror("set_vector is only for vectors.");
	int vectorindex = ptTypeSpecificIndex[genindex];
	int veclength = ptIntSizeOfVectorParameters[vectorindex];
	if (i < 0 || i >= veclength) srerror("out of range in get_vector_val.");
	(ptptValuesOfVectorParameters[vectorindex])[i] = v;
};
