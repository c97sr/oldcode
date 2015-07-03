# First run to see if the implications of the US R0, the Oz case mix and the QLD data

source("age_dynamics_setup.R")

mcmcH1N1v1(	"fitted_params_inc_r.csv",
			"run_1.csv",
			FourAgeParams(),
			incData,
			PdmSolveLikeV2,
			case_mix=c(7,81,118,15),			# Australian data
			# case_mix=c(228,1194,840,25),		# Early US data
			no_samples=10,
			samp_freq=1,
			report_block=2)
	