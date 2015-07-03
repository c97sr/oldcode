#include "sim.h"

/*
 * Current notes for development
 * - the main task is to construct a system with arbitrary spatial scale
 * - the smallest spatial scale is the landscan resolution
 * - the system will run at that scale, so simulation will be on the landscan units
 * - however the actual distances used will be determined by a less highly resolved scale
 * - output variables of interest will be timeseries of incidence across populations
 * - as bits start to work, make sure that the methdos sections are kept uptodate
 *
 * Next
 * - Use the R code postproc.R to debug the foi calculation for this simulation
 * -
 *
 *
 *
 */

// Use the standard name space
using namespace std;

// Start the main loop. No cammand line arguments at the moment
int main() {

	// Declare some useful "global" constants
	double epsilon = 1e-30;

	// Read the parameters, population sizes and their connectivities
	Parameters p("parameters.csv");
	string str_fn_pop_dens = p.get_string_val("fn_ascii_pop");

	// Declare and initialize (if needed) variables required in the main scope
	double dbl_dt = 	static_cast<double>(p.get_scalar_val("dt"));
	double dbl_dtm = 	static_cast<double>(p.get_scalar_val("dtm"));
	double dbl_tf = 	static_cast<double>(p.get_scalar_val("tf"));
	double dbl_D_I = 	static_cast<double>(p.get_scalar_val("D_I"));
	double dbl_beta = 	static_cast<double>(p.get_scalar_val("R_0"))/dbl_D_I;
	double dbl_D_R = 	static_cast<double>(p.get_scalar_val("D_R"));
	double dbl_beta_s = static_cast<double>(p.get_scalar_val("beta_s"));

	// Construct the matrix to be used for output for a single realization
	// Potential problem here if dt_m is not a multiple of dt
	int intNoSimTimes = 	static_cast<int>(dbl_tf/dbl_dt)+1;
	int intNoMeasureTimes = static_cast<int>(dbl_tf / dbl_dtm+1);
	int int_nr = 			static_cast<int>(p.get_scalar_val("nr"));

	// Read in the ascii grid file header info
	AsciiPop gzpop(str_fn_pop_dens,1);

	// Derive the required parameters from the AsciiPop
	int noBasicX = gzpop.GetNoX();
	int noBasicY = gzpop.GetNoY();
	int noBasicA = 1;

	// Construct any auxiliary data structures required for the simulation
	gsl_matrix * matCumInc = gsl_matrix_alloc(intNoMeasureTimes , noBasicX * noBasicY);	// The main output data structure for incidence
	ConnectMatrix connect(noBasicA,noBasicX,noBasicY);											// The connectivity matrix between two different populations

	// Declare and initialize the state variables and "force of" variables
	StateVariable sus(gzpop);
	StateVariable pop(gzpop);
	StateVariable inf(gzpop);
	StateVariable rec(gzpop);
	StateVariable cuminc(gzpop);
	StateVariable foinf(gzpop);
	StateVariable forec(gzpop);
	StateVariable fobsa(gzpop);
	StateVariable evInf(gzpop);
	StateVariable evRec(gzpop);
	StateVariable evBsa(gzpop);

	// Set up the machinery for the random number generation
	const gsl_rng_type * T;
	gsl_rng * r;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc (T);

	// Declare some variables required in the main loops
	double t,tmpval;

	// Start the main replication loop
	for (int i=0;i<int_nr;i++) {

		// Initialize the state variables
		sus.ResetOneAge(gzpop);
		inf.setall(0);
		rec.setall(0);
		cuminc.setall(0);

		// Start the main time loop
		int intNextMeasureIndex = 0;
		for (int j=0;j<intNoSimTimes;j++) {

			// Initialize variables required for the timestep loop
			t = static_cast<double>(j)*dbl_dt;
			foinf.setall(0);
		 	forec.setall(0);
		 	fobsa.setall(0);

			// Log progress to the screen
			cout << "Time " << t << " of " << dbl_tf << " in realization " << i+1 << " of " << int_nr << ".       \r";

			// Test to see if the next measure time has been passed and
			// record output values if it has
			if (t - intNextMeasureIndex*dbl_dtm + epsilon > 0) {

				// Put cumulative numbers of infections into the reporting matrix
				for (int k=0;k<noBasicX;k++) {
					for (int l=0;l<noBasicY;l++) {
						tmpval = cuminc.get(0,k,l);
						gsl_matrix_set(matCumInc,intNextMeasureIndex,k*noBasicY+l,tmpval);
					}
				}

				// Increment the measure time index
				intNextMeasureIndex++;

			// Close the bracket for measure times
			}

			// Calculate the hazards of infection, recovery and becoming susceptible
			// Just not working well at all
			calc_foi(foinf,sus,inf,rec,pop,connect,dbl_beta,dbl_beta_s);
			calc_genhaz(forec,inf,1.0 / dbl_D_I);
			calc_genhaz(fobsa,rec,1.0 / dbl_D_R);

			// Generate the events
			gen_genev(evInf, sus, foinf,r,dbl_dt);
			gen_genev(evRec, inf, forec,r,dbl_dt);
			gen_genev(evBsa, rec, fobsa,r,dbl_dt);

			// Update the state variables
			sus.add(evBsa,1);
			sus.add(evInf,-1);
			inf.add(evInf,1);
			inf.add(evRec,-1);
			rec.add(evRec,1);
			rec.add(evBsa,-1);
			cuminc.add(evInf,1);

		}

		// Write the file for a single realization, file naming and format are designed to be read into R as a 4D array
		write_realization_to_file(matCumInc,"testout",i,noBasicX,noBasicY);

	// End the per realization loop
	}

	// End the logging line
	cout << "\n";

	// Free up anything explicitly allocated in the main scope
	gsl_matrix_free(matCumInc);
	gsl_rng_free (r);

	// Return 0 to indicate successful completion. Change here to check compilation string.
	cout << "Program completed normally." << endl;
	return 0;

// End main()
}


