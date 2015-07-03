FourAgeParams <- function(){
	
	c(		N 			= 1,
			seed 		= 1000/7000000,
			r 			= 0.26,
			Tg 			= 2.6,
			phi_1		= 1.0,
			phi_2		= 1.0,
			phi_3		= 1.0,
			phi_4		= 1.0,
			p_R 		= 0.86,
#			p_H_base 	= 0.05,
#			p_H_base_1 	= 0.8233,
#			p_H_base_2 	= 0.1992,
#			p_H_base_3 	= 1.0,
#			p_H_base_4 	= 0.2143,
			p_H_base 	= 1,
			p_H_base_1 	= 1,
			p_H_base_2 	= 1,
			p_H_base_3 	= 1,
			p_H_base_4 	= 1,
# XX
			p_I 		= 0.13,
			p_I_1 		= 0.13,
			p_I_2 		= 0.13,
			p_I_3 		= 0.13,
			p_I_4 		= 0.13,
			gamma_v 	= 1/7.0,							# NEJM ANZ piece 10.1056/NEJMoa0908481
			# gamma_h 	= 1/7.0,							# Feagan
			gamma_h 	= 1/3.0,							# Feagan
			t0 			= as.numeric(as.date("21May2009")),
			tf			= as.numeric(as.date("31Dec2009")),
			# From http://www.census.gov/popest/national/asrh/NC-EST2008-sa.html
			# From http://www.oesr.qld.gov.au/queensland-by-theme/demography/population/tables/pop-proj-medium/proj-pop-age-qld/index.shtml
			p_1			=  268090,	
			p_2			=  794098,								# Weighted probability factors
			p_3			= 2535833,
			p_4 		=  493522,
			# p_1			=  	5,	
			# p_2			=  	12,								# Weighted probability factors
			# p_3			= 	50,
			# p_4 		=  	10,
			mixmatindex	= 2,									# Select
			fitsus = 0,
			t_st 		= as.numeric(as.date("30June2009")),
			delta_sm	= 1.0
	)	
}
