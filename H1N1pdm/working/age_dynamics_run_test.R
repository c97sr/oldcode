# First run to see if the implications of the US R0, the Oz case mix and the QLD data

source("age_dynamics_setup.R")

set.seed(1234)

mcmcH1N1v1(	"fitted_params_inc_r.csv",
			"test_mcmc.csv",
			FourAgeParams(),
			PdmSolveLikeV2,
			no_samples=10,
			samp_freq=1,
			report_block=10,
			casemix=c(228,1194,840,25),
			obs=incData,
			ngm=MakeNgmFour
)		

# read.csv("test_mcmc.csv")