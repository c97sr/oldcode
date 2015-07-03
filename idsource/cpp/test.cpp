#include "sim.h"

int main() {

	// Current objectives
	// - Define and load up all the inputs as they are defined in the methods section
	// - Define and export all the outputs required for the analysis
	// - Sort all the house keeping

	// Read in a list of populations with age breakdowns
	Parameters parameters("parameters.csv");
	
	
	
	cout << parameters.get_string_val("fnPopdef") << endl;
	string * ptTmpStr = parameters.get_string_pt("fnPopdef");
	cout << *ptTmpStr << endl;
	parameters.set_string("fnPopdef","alibaba");
	cout << parameters.get_string_val("fnPopdef") << endl;
	cout << parameters.get_scalar_val("alpha") << endl;
	double * ptTmpDbl = parameters.get_scalar_pt("alpha");
	cout << *ptTmpDbl << endl;
	parameters.set_scalar("alpha",34.5);
	cout << *ptTmpDbl << endl;
	cout << parameters.get_vector_val("t_m",1) << endl;
	parameters.set_vector("t_m",1,12);
	cout << *(parameters.get_vector_pt("t_m")+1) << endl;
	parameters.set_string("fnPopdef","m_test.csv");
	gsl_matrix * mat = gsl_matrix_read_csv(parameters.get_string_val("fnPopdef"));
	gsl_matrix_set(mat,1,1,100);
	gsl_matrix_write_csv(mat,"test_out.csv");
	gsl_matrix_free(mat);

	return 0;

}


