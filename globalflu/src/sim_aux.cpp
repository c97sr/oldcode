#include"sim.h"

// Utility error function
// Takes a
void srerror(const string etext) {
	cerr << "Run-time error..." << endl;
	cerr << etext << endl;
	cerr << "...now exiting." << endl;
	exit(1);
};

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

gsl_matrix * gsl_matrix_read_csv(string filename){

	// Read in the values from filename and return a gsl matrix
	// Assumes some fairly strict formatting of the text file

	// Declare some required local variables
	string filecontents,tmp1,tmp2;
	int norows,nocols,rowlength;

	// Take the contents of the file
	ifstream ifs(filename.c_str());
	if (ifs.fail()) srerror("Can't open the file in constructor for matrix read: "+filename);
	filecontents = "";
	while (ifs) {
		ifs >> tmp1;
		filecontents = filecontents + tmp1 + "\n";
	}
	ifs.close();

	// Cycle through to measure the size of the matrix and make sure its rectangular
	int currentrow = 0;
	nocols = 0;
	stringstream iss1(filecontents);
	while (getline(iss1,tmp1)) {
		if (tmp1 != "") {
			stringstream ifs2(tmp1);
			nocols = 0;
			while (getline(ifs2,tmp2,',')) {
				nocols++;
			}
			if (currentrow==0)rowlength = nocols;
			else if (nocols != rowlength) srerror("matrix seems not to be rectangular");
			currentrow++;
		}
	}
	norows = currentrow - 1;

	// Declare the gsl matrix
	gsl_matrix * m = gsl_matrix_alloc (norows, nocols);

	// Read in the values of the matrix
	double dbltmp;
	stringstream iss2(filecontents);
	for (int i=0;i<norows;i++) {
		for (int j=0;j<nocols;j++) {
			iss2 >> dbltmp;
			gsl_matrix_set(m,i,j,dbltmp);
		}
	}

	// return the pointer to the matrix
	return m;

	// end function
};

int gsl_matrix_write_csv(gsl_matrix * m, string filename) {

	// Write the contents of m to file filename
	ofstream s(filename.c_str());
	s.precision(6);
	if (s.fail()) return 0;
	for (unsigned int i=0;i<m->size1;i++) {
		unsigned int j=0;
		for (j=0;j < (m->size2-1);j++)
			s << gsl_matrix_get(m,i,j) << ",";
		s << gsl_matrix_get(m,i,j) << "\n";
	}

	// Close the stream and return a success signal
	s.close();
	return 1;

};

AsciiPop::AsciiPop(string filename, bool report) {

	static int nullData;
	ifstream ifs;
	double tmp,total=0;
	string junk;
	vector<double>::iterator ptVecDbl;
	ifs.open(filename.c_str());
	if (ifs.fail()) srerror("Problem opening density file");

	ifs >> junk >> nox;
	ifs >> junk >> noy;
	ifs >> junk >> minx;
	ifs >> junk >> miny;
	ifs >> junk >> stepx;
	ifs >> junk >> nullData;
	minx = minx / 57.2957795131;
	miny = miny / 57.2957795131;
	stepx = stepx / 57.2957795131;
	stepy=stepx;


	vals = new double[nox*noy]; for (int i=0;i<nox*noy;++i) vals[i]=0;

	maxx=minx+(nox-1)*stepx;
	maxy=miny+(noy-1)*stepy;

	for (int i=noy-1;i>=0;--i) {
		for (int j=0;j<nox;++j) {
			ifs >> tmp;
			if (tmp == nullData) tmp = 0;
			// vals[j*noy+(noy-i-1)]=tmp;
			vals[j*noy+i]=tmp;
			total += tmp;
		}
	}

	ifs.close();
	CalcMaxVal();
	if (report) cerr << "Ascii pop constructed.  Total: " << total << "\n";
};

string AsciiPop::Table() {
	ostringstream oss;
	oss << "0\t";
	for (int i=0;i<nox;++i) oss << minx+1.0*i*stepx << "\t";
	oss << "\n";
	for (int j=0;j<noy;++j) {
		oss << 1.0*j*stepy+miny << "\t";
		for (int i=0;i<nox;++i) oss << vals[i*noy+j] << "\t";
		oss << "\n";
	}
	return oss.str();
};

double AsciiPop::Value(double x, double y) {
	double rtnval;
	static int xcoord,ycoord;
	if (x < minx || x > maxx || y < miny || y > maxy) srerror("coords out of range in DensityField::Value(double x, double y)");
	xcoord = static_cast<int>((x-minx)/stepx);
	ycoord = static_cast<int>((y-miny)/stepy);
	rtnval = (vals[xcoord*noy+ycoord]+vals[(xcoord+1)*noy+ycoord]+vals[(xcoord+1)*noy+ycoord+1]+vals[xcoord*noy+ycoord+1])/4;
	return rtnval;
};

double AsciiPop::Value(int x, int y) {
	double rtnval;
	static int xcoord,ycoord;
	if (x < 0 || x >= nox || y < 0 || y >= noy) srerror("coords out of range in DensityField::Value(int x, int y)");
	xcoord = x;
	ycoord = y;
	rtnval = vals[xcoord*noy+ycoord];
	return rtnval;
};

void AsciiPop::CalcMaxVal() {
	double* pt = vals;
	double rtnval = -9999;
	while (pt != vals+nox*noy) {
		if (*pt > rtnval) rtnval = *pt;
		pt++;
	}
	if (rtnval==0) srerror("Entire density field is zero.");
	maxval = rtnval;
};

ConnectMatrix::ConnectMatrix(int noa_in, int nox_in, int noy_in) {

	// Declare the 4 d array
	nox = nox_in;
	noy = noy_in;
	noa = noa_in;
	cons = gsl_matrix_alloc(noa*nox*noy,noa*nox*noy);
	if (cons==NULL) srerror("Problem allocating memory in ConnectMatrix::ConnectMatrix(int nox_in, int noy_in)");

	// Initialize all the values to 1
	for (int i=0;i<noa;i++)
		for (int j=0;j<nox;j++)
			for (int k=0;k<noy;k++)
				for (int l=0;l<noa;l++)
					for (int m=0;m<nox;m++)
						for (int n=0;n<noy;n++)
								gsl_matrix_set(cons,i*nox*noy+j*noy+k,l*nox*noy+m*noy+n,1);

};

ConnectMatrix::~ConnectMatrix() {

	// Free any allocated memory for the object
	gsl_matrix_free(cons);

};

void write_realization_to_file(gsl_matrix * m, string filestem, int real_number, int nox, int noy) {

	// Function to take cumulative incidence matrix, the number of x squares and the number of y squares and the time step and
	// output a csv file that r can read to reconstruct an xyz array for the data

	// Set up the correct filename
	stringstream ss;
	ss << real_number;
	string filename = filestem + ss.str() + ".csv";

	// Write the contents of the matrix and some housekeeping items to file filename
	ofstream s(filename.c_str());
	s.precision(6);
	if (s.fail()) srerror("Unable to write in write_realization_to_file");
	for (unsigned int i=0;i< m->size1;i++) {
		s << i << "," << nox << "," << noy ;
		for (unsigned int j=0;j < m->size2;j++) s << ","  << gsl_matrix_get(m,i,j);
		s << "\n";
	}

	s.close();
};

StateVariable::StateVariable(AsciiPop &ap) {

	nox = ap.GetNoX();
	noy = ap.GetNoY();
	noa = 1;

	vals = gsl_vector_alloc(noa*nox*noy);
	if (vals==NULL) srerror("Problem allocating memory in ConnectMatrix::ConnectMatrix(int nox_in, int noy_in)");
	for (int i=0;i<noa;i++)
		for (int j=0;j<nox;j++)
			for (int k=0;k<noy;k++)
				gsl_vector_set(vals,i*nox*noy+j*noy+k,ap.Value(j,k));

};

void StateVariable::ResetOneAge(AsciiPop &ap) {

	for (int j=0;j<nox;j++)
		for (int k=0;k<noy;k++)
			gsl_vector_set(vals,j*noy+k,ap.Value(j,k));

};

StateVariable::~StateVariable(){

	gsl_vector_free(vals);

};

void calc_foi(	StateVariable &foi_in, StateVariable &sus_in, StateVariable &inf_in, StateVariable &rec_in, StateVariable &pop_in,
				ConnectMatrix &cm_in, double beta, double beta_s) {

	// Declare the variables
	static double epsilon = 1e-30;
	double tmpnosus, tmpnoinf, tmpfoi, popsource, poptarget;
	int nox, noy, noa;
	nox = foi_in.getnox();
	noy = foi_in.getnoy();
	noa = foi_in.getnoa();

	// Set up the outer part of the six way loop, for the susceptible sweep
	for (int i=0;i<noa;i++) {
		for (int j=0;j<nox;j++) {
			for (int k=0;k<noy;k++) {

				// Check to see if the inner infectious loop is required
				tmpnosus = sus_in.get(i,j,k);
				if (tmpnosus - epsilon > 0) {

					// Initialize the force of infection felt by the susceptibles
					tmpfoi = beta_s;

					// If there are susceptibles, need to check out the infectious
					for (int l=0;l<noa;l++) {
						for (int m=0;m<nox;m++) {
							for (int n=0;n<noy;n++) {

								// Assign the number of infectious
								tmpnoinf = inf_in.get(l,m,n);

								// Up to here. Think this below is just bollocks!
								// Need to make the force of infection reflect that m is a proportion of time
								// This will require a static effective population size for each cell square that incorporates travel

								// Only add to the foi if there are infectious individuals present in the target square
								if (tmpnoinf - epsilon > 0) {

									// Access the required variables
									popsource = pop_in.get(i,j,k);
									poptarget = pop_in.get(l,m,n);

									// Add to the required increment to the force of infection
									tmpfoi += beta / 2.0 * (	cm_in.get(i,j,k,l,m,n)*tmpnoinf/poptarget +
																cm_in.get(l,m,n,i,j,k)*tmpnoinf/popsource );

								}
							}
						}
					}

					// Set the force of infection for the source population
					foi_in.set(i,j,k,tmpfoi);

				// Else statement for the test of presence of susceptibles
				} else {

					// If there are no susceptible individuals, add in an error value to the foi
					foi_in.set(i,j,k,-1);

				}
			}
		}
	}

};

// Takes a state variable s and a rate p and generates a matrix of hazards that the event occurs
void calc_genhaz(StateVariable &h, StateVariable &s, double p) {

	// Declare the variables
	static double tmpno;
	static int nox, noy, noa;
	nox = s.getnox();
	noy = s.getnoy();
	noa = s.getnoa();

	// Set up the 3D loops
	for (int i=0;i<noa;i++) {
		for (int j=0;j<nox;j++) {
			for (int k=0;k<noy;k++) {

				// Extract the population value and multiply the hazard
				tmpno = s.get(i,j,k)*p;
				h.set(i,j,k,tmpno);

			}
		}
	}

};

// Takes a state variable and a known hazard per sub pop and a time step and
// fills e with a random set of binomial draws
void gen_genev(StateVariable &e, StateVariable &s, StateVariable &h, gsl_rng * r, double dt) {

	// Declare the aux variables
	static double epsilon = 1e-20;
	static double tmpno, tmphaz, tmpprob;
	static unsigned int tmppop;
	static int nox, noy, noa;
	nox = s.getnox();
	noy = s.getnoy();
	noa = s.getnoa();

	// Setup the vectors for the random number call
	static size_t K = 2;
	static double p[2];
	static unsigned int n[2];

	// Set up the 3D loops
	for (int i=0;i<noa;i++) {
		for (int j=0;j<nox;j++) {
			for (int k=0;k<noy;k++) {

				// Check for nonzero hazard and non zero population
				tmppop = static_cast<unsigned int>(s.get(i,j,k));
				tmphaz = h.get(i,j,k);
				if (tmppop - epsilon > 0 && tmphaz - epsilon > 0) {

					// Calculate the probability
					tmpprob = 1-exp(-tmphaz*dt);
					p[0] = tmpprob;
					p[1] = 1- tmpprob;
					gsl_ran_multinomial(r,K,tmppop,p,n);
					tmpno = static_cast<double>(n[0]);
					e.set(i,j,k,tmpno);

				} else {

					// If either the hazard of the state variable are zero then there can be no events
					e.set(i,j,k,0);

				}

			}
		}
	}

};

// Adds all the members of one state variabel to another
void StateVariable::add(StateVariable &sv, double factor) {

	// Define some variables for the function
	double argnoa,argnox,argnoy,tmparg,tmpthis,val;

	// Record properties of the argument StateVariable
	argnox = sv.getnox();
	argnoy = sv.getnoy();
	argnoa = sv.getnoa();

	// Check for compatibility
	if (argnox != nox) srerror("problem in StateVariable::add");
	if (argnoy != noy) srerror("problem in StateVariable::add");
	if (argnoa != noa) srerror("problem in StateVariable::add");

	// Set up the 3D loops
	for (int i=0;i<noa;i++) {
		for (int j=0;j<nox;j++) {
			for (int k=0;k<noy;k++) {

				// Extract the population value and multiply the hazard
				tmparg = sv.get(i,j,k);
				tmpthis = get(i,j,k);
				val = tmparg+tmpthis;
				this->set(i,j,k,val);

			}
		}
	}

};

// Write a state variable to the standard error
void StateVariable::PrintCerr(int age) {
	for (int i=0;i<nox;i++) {
		for (int j=0;j<noy;j++) {
			cerr << this->get(age,i,j) << "\t";
		}
		cerr << "\n";
	}
	cerr << "\n";q
}
