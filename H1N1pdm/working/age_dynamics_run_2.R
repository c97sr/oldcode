# First run to see if the implications of the US R0, the Oz case mix and the QLD data

source("age_dynamics_setup.R")

mcmcH1N1v1(	"fitted_params_inc_r.csv",
			"run_2.csv",
			FourAgeParams(),
			PdmSolveLikeV2,
			obs=incData,
			# case_mix=c(7,81,118,15),			# Australian data
			casemix=c(228,1194,840,25),		# Early US data
			no_samples=8000,
			samp_freq=10,
			report_block=5,
			ngm=MakeNgmFour)	