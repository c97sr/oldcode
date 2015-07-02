#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<iostream>
#include<vector>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include<cmath>

using namespace std;

// Little utility function for errors
void srerror(const string etext);

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

class AsciiPop {
private:
	double* vals;
	double minx,miny,maxx,maxy,stepx,stepy,maxval;
	int nox,noy; // number of grid points
	void CalcMaxVal();
public:
	// nox and noy are numbers of intervals in the constructors, not number of points
	AsciiPop() {srerror("No default constructor for AsciiPop");};
	// DensityField(double d, int nox_in, int noy_in, double minx_in, double maxx_in, double miny_in, double maxy_in);
	// DensityField(double f(double,double), int nox_in, int noy_in, double minx_in, double maxx_in, double miny_in, double maxy_in);
	// DensityField(string filename, double stepx_in, double stepy_in, double minx_in, double maxx_in, double miny_in, double maxy_in);
	AsciiPop(string filename, bool report);
	// ~DensityField(){delete [] vals;};
	~AsciiPop(){};
	double Value(double x, double y);
	double Value(int x, int y);
	string Table();
	inline double GetMinX(){return minx;};
	inline double GetMaxX(){return maxx;};
	inline double GetMinY(){return miny;};
	inline double GetMaxY(){return maxy;};
	inline double GetMaxVal(){return maxval;};
	inline int GetNoX(){return nox;};
	inline int GetNoY(){return noy;};
};

class ConnectMatrix {
private:
	int nox,noy,noa;
	gsl_matrix * cons;
public:
	ConnectMatrix() {srerror("No default constructor for ConnectMatrix");};
	ConnectMatrix(int nox_in, int noy_in, int noa_in);	// Constructs a completely flat, mass action kernel
	~ConnectMatrix();
	inline double get(int a1, int x1, int y1, int a2, int x2, int y2) {return gsl_matrix_get(cons,a1*nox*noy+x1*noy+y1,a2*nox*noy+x2*noy+y2);};
};

// An structure for arbitrary age classes in the metapopulations
class StateVariable {
private:
	int nox,noy,noa;
	gsl_vector * vals;
public:
	StateVariable() {srerror("No default constructor for ConnectMatrix");};
	StateVariable(AsciiPop &ap);	// Constructs a single age class population from a single ascii pop
	~StateVariable();
	void add(StateVariable &sv, double factor);
	void ResetOneAge(AsciiPop &ap);	// Sets the values equal to the asciipop values
	void PrintCerr(int age); // Sends the contents of the State Variable to the standard error
	inline double get(int a, int x, int y) {return gsl_vector_get(vals,a*nox*noy+x*noy+y);};
	inline double getallage(int x, int y) {static double rtn = 0; for (int i=0;i<noa;i++) rtn += get(i,x,y); return rtn;}
	inline void set(int a, int x, int y, double v) {gsl_vector_set(vals,a*nox*noy+x*noy+y,v);};
	inline void setallage(int a, double v) {for (int i=0;i<nox;i++) for (int j=0;j<noy;j++) set(a,i,j,v);};
	inline void setall(double v) {for (int a=0;a<noa;a++) setallage(a,v);};
	inline int getnox() {return nox;};
	inline int getnoy() {return noy;};
	inline int getnoa() {return noa;};
	inline double sum() {static double rtn = 0; for (int i=0;i<noa;i++) for (int j=0;j<nox;j++) for (int k=0;k<noy;k++) rtn += get(i,j,k); return rtn;}
};

// Functions that don't fit neatly in a single class
gsl_matrix * gsl_matrix_read_csv(string filename);
int gsl_matrix_write_csv(gsl_matrix * m, string filename);
void write_realization_to_file(gsl_matrix * cuminc, string filestem, int real_number, int nox, int noy);
void calc_foi(StateVariable &foi_in, StateVariable &sus_in, StateVariable &inf_in, StateVariable &rec_in, StateVariable &pop_in,
		ConnectMatrix &cm_in, double beta, double beta_s);
void calc_genhaz(StateVariable &h, StateVariable &s, double p);
void gen_genev(StateVariable &e, StateVariable &s, StateVariable &h, gsl_rng * r, double dt);
