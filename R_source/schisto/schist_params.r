parvector =  
    c(  
        alpha_S = 1/((11.18+6.21)/2/52),    # Anderson and May p 558, Pesigan 1958
        sigma_S = 1/((70+53)/2/365),		# Anderson and May p 563, Pesigan 1958, Hairston 1965a
        beta_SM = 1,            			# Auxilliary
        beta_MS = 1,            			# Auxilliary
        gamma_H = 0,          				# To be estimated
        gamma_L = 0,          				# To be estimated 
        epsilon_D = 0.0001,     			# To be estimated
        epsilon_C = 0.0001,     			# To be estimated
        epsilon_P = 0.0001,     			# To be estimated
        epsilon_W = 0.001,      			# To be estimated
        epsilon_R = 0.001,      			# To be estimated
        epsilon_L = 1,          			# Set to unity wihtout loss of generality.
        epsilon_H = 0.01,       			# To be estimated
        omega_D = 0.0001,       			# To be estimated
        omega_C = 0.0001,       			# To be estimated
        omega_P = 0.0001,       			# To be estimated
        omega_W = 0.001,        			# To be estimated
        omega_R = 0.001,        			# To be estimated
        omega_L = 1,            			# Set to unity wihtout loss of generality.
        omega_H = 0.01,         			# To be estimated
        gamma_D = 1/(20.8/12),  			# Fitting exponential lifespan to ages in study
        gamma_C = 1/(19.3/12),  			# Fitting exponential lifespan to ages in study
        gamma_P = 1/(5.9/12),   			# Fitting exponential lifespan to ages in study
        gamma_W = 1/(96/12),    			# Distribution more like type 2 survival, 95% of WB younger than this age
        gamma_R = 1/(12/12),    			# Approximate guess
        N_S = 1,	             			# Should be N_S=1
        N_H = 1000,             			# Entered from study data
        N_P = 1000,             			# Entered from study data
        N_C = 1000,             			# Entered from study data
        N_D = 1000,             			# Entered from study data
        N_W = 1000,             			# Entered from study data
        N_R = 1000              			# Entered from study data
)
