# Run the four-class transmission model with the MLE parameters from the MCMC runs
baseHK4 <- 
	c(		N 			= 7000000,
			seed 		= 1,
			r 			= 0.074169115,
			Tg 			= 2.6,
			phi_1		= 6.16098487,
			phi_2		= 7.415206875,
			phi_3		= 1.0,
			phi_4		= 0.274544092,
			p_R 		= 1.0,
			p_H_base 	= 0.047216386,
			p_H_base_1 	= 0.816666667,
			p_H_base_2 	= 0.203333333,
			p_H_base_3 	= 1.0,
			p_H_base_4 	= 0.613333333,
			p_I 		= 1.0,
			p_I_1 		= 1.0,
			p_I_2 		= 1.0,
			p_I_3 		= 1.0,
			p_I_4 		= 1.0,
			gamma_v 	= 1/7.0,							# NEJM ANZ piece 10.1056/NEJMoa0908481
			gamma_h 	= 1/3.0,							# Feagan
			t0 			= 18013.20894,
			tf			= as.date("31Mar2010"),
			p_1			= 0.05,	
			p_2			= 0.15,								# Weighted probability factors
			p_3			= 0.4,
			p_4 		= 0.4,
			mixmatindex	= 2,								# Select
			fitsus		= 0,
			t_st		= as.date("01Sep2009"),
			delta_sm	= 1									# Reduction in mixing between school age children during school holidays
	)	
