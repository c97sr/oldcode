## Things to be changed
no_samples = 1
hc_lb = c(1.1,0.2,0.2,0.5,0.2)
hc_ub = c(4,0.5,0.3,0.75,0.5)
single_set = c(1.8,0.35,0.25,0.666,0.5)
scale_factors = c(1,1,1,1,1/100000)
hc_log = c(FALSE,FALSE,FALSE,FALSE,FALSE)
opt_lower_bounds = c(0.00001,0.00001,0.00001,0.0001,0.00000000001)
lower_bounds = c(0.01,0.01,0.01,0.01,0.00001)
upper_bounds = c(10,5,5,5,0.001)
no_iterations = 50
outfile="D:\\share\\files\\projects\\flu2005\\city_contact\\param_sets\\R_script_HC.out"

# R0 = average number of secondary cases
# Theta = Proportion of transmission which is presymptomatic or asymptomatic
# Omega = Proportion of theta that occurs as asymptomatic transmission
# Alpha = Proportion of secondary cases that occur outside the home
# Gamma = Proportion of transmission outside the home due to community transmission
# v[1] = basic network beta
# v[2] = rel hazard presymptomatic to post symptomatic
# v[3] = relative hazard asymtpomatic to presymptomatic
# v[4] = relative hazard in the peer group compared to the home
# v[5] = relative hazard in the community versus the peer group 

p_H_basehs = array(c(0.0599,0.1608,0.2401,0.3188,0.1519,0.0513,0.0118,0.0035,0.0010,0.0010),c(10))
p_ad = 5027968 / (5027968 + 1213322)

max_housemates = dim(p_H_basehs)[1]		## max number of other members of household
mu_P = 0.5						## average duration of presymptomatic stage		
mu_S = 3						## Average duration of symptomatic stage
p_S = 0.68						## Proportion of infections which are symptomatic
n_WP_raw = 30					## Expected number of workplace collegues (uniform 
n_SC_raw = 40					## Expected number of school class mates
delta_SP_A = 0.1					## Probabilility of attending work when symptomatic
delta_SP_C = 0.1					## Probability of attending school when symptomatic

n_WP = n_WP_raw+1-1
var_SC = 1/12*(45-35)^2
n_SC = (var_SC+n_SC_raw^2)/n_SC_raw-1 

## Make a latin hypercube
unit_to_log <- function(v,lb,ub) {
	rtnval = lb*10^((log10(ub)-log10(lb))*v)
	rtnval
}

unit_to_linear <- function(v,lb,ub) {
	rtnval = lb + (ub-lb)*v
	rtnval
}

hyper_vector <- function(n,lb,ub,log) {
	wk = array(runif(n,0,1),c(n))
	for (i in 1:n) {
		wk[i] = (i-1)/n+wk[i]/n
		if (log) wk[i] = unit_to_log(wk[i],lb,ub) 
		else wk[i] = unit_to_linear(wk[i],lb,ub)
	}
	ranrank = array(runif(n,0,1),c(n))
	index = order(ranrank)
	rtnval = array(0,c(n))
	for (i in 1:n) rtnval[i]=wk[index[i]]
	rtnval
}

phi_P <- function(hazard,meantime) {
	rtnval <- 1-exp(-hazard*meantime)
	rtnval
}

phi_S <- function(haz,mt,alph) {
	b = mt/alph
	phi_S_rtnval <- 1 - 1 /(1+haz*b)^alph
	phi_S_rtnval
}

R_PH <- function(haz) {
	rtnval = 0
	for (i in 1:max_housemates) {
		rtnval = rtnval + p_H_basehs[i] * (i-1) * phi_P(haz/i,mu_P)  
	}
	rtnval
}

R_SH <- function(haz,prehaz) {
	rtnval = 0
	for (i in 1:max_housemates) {
		rtnval = rtnval + p_S * p_H_basehs[i] * (1 - phi_P(prehaz/i,mu_P)) * (i - 1) * phi_S(haz/i,mu_S,1)   
	}
	rtnval
}

R_AH <- function(haz,prehaz) {
	rtnval = 0
	for (i in 1:max_housemates) {
		rtnval = rtnval + (1 - p_S) * p_H_basehs[i] * (1 - phi_P(prehaz/i,mu_P)) * (i-1) * phi_S(haz/i,mu_S,1)   
	}
	rtnval
}

R_PP <- function(haz) {
	rtnval = 0
	rtnval = rtnval + n_WP * phi_P(haz,mu_P)  
	rtnval
}

R_SP <- function(haz,prehaz) {
	rtnval = 0
	rtnval = rtnval + delta_SP_A * p_ad * p_S * (1 - phi_P(prehaz,mu_P)) * n_WP * phi_S(haz,mu_S,1)   
	rtnval = rtnval + delta_SP_C * (1 - p_ad ) * p_S * (1 - phi_P(prehaz,mu_P)) * n_SC * phi_S(haz,mu_S,1)   
	rtnval
}

R_AP <- function(haz,prehaz) {
	rtnval = 0
	rtnval = rtnval + (1 - p_S) * (1 - phi_P(prehaz,mu_P)) * n_WP * phi_S(haz,mu_S,1)   
	rtnval
}


## r_PC = N*(1-PMGF(presymHazard*infHazardRandom*beta));
## r_AC = (1-pSymptoms)*N*PMGF(presymHazard*infHazardRandom*beta)*(1-IMGF(asymHazard*infHazardRandom*beta));
## r_SC = pSymptoms*N*PMGF(presymHazard*infHazardRandom*beta)*(1-IMGF(infHazardRandom*beta));

reparam <- function(v) {
	popsize=6800000
	# R0 = average number of secondary cases
	# Theta = Proportion of transmission which is presymptomatic or asymptomatic
	# Omega = Proportion of transmission those never symptomatic
	# Alpha = Proportion of secondary cases that occur outside the home
	# Gamma = Proportion of transmission outside the home due to community transmission
	# v[1] = basic network beta
	# v[2] = rel hazard presymptomatic to post symptomatic
	# v[3] = relative hazard asymtpomatic to symptomatic
	# v[4] = relative hazard in the peer group compared to the home
	# v[5] = relative hazard in the community versus the peer group 
	if (min(v) < 0) {rtnval = 1e100} 
	else {
		r_PH = R_PH(v[2]*v[1])
		r_SH = R_SH(v[1],v[2]*v[1])
		# r_AH = R_AH(v[1]*v[2]*v[3],v[2]*v[1])
		r_AH = R_AH(v[1]*v[3],v[2]*v[1])
		r_PP = R_PP(v[1]*v[2]*v[4])
		r_SP = R_SP(v[1]*v[4],v[1]*v[2]*v[4])
		r_AP = R_AP(v[1]*v[3]*v[4],v[2]*v[1]*v[4])
		r_PC = popsize*v[1]*v[2]*v[5]*mu_P
		r_SC = popsize*p_S*v[1]*v[5]*mu_S
		r_AC = popsize*(1-p_S)* v[1]*v[3]*v[5]*mu_S
		# r_AP = R_AP(v[1]*v[2]*v[3]*v[4],v[2]*v[1]*v[4])
		## r_AC = popsize*(1-p_S)* v[1]*v[2]*v[3]*v[4]*v[5]*mu_S	
		## r_PC = popsize*phi_P(v[1]*v[4]*v[5],mu_P)
		## r_SC = popsize*p_S * phi_S(v[1]*v[4]*v[5],mu_S,1)*(1-phi_P(v[1]*v[2]*v[4]*v[5],mu_P))
		## r_AC = popsize*(1-p_S)*phi_S(v[1]*v[2]*v[3]*v[4]*v[5],mu_S,1)*(1-phi_P(v[1]*v[2]*v[4]*v[5],mu_P))	
		R0 = r_PH + r_SH + r_AH + r_PP + r_SP + r_AP + r_PC + r_AC + r_SC
		Theta = (r_PH + r_AH + r_PP + r_AP + r_PC + r_AC) / R0
		# Omega = (r_AH + r_AP + r_AC) / (r_AH + r_PH + r_AP + r_PP + r_AC + r_PC)
		Omega = ((r_AH + r_AP + r_AC) + (1-p_S)* (r_PH + r_PP + r_PC)) / R0
		Alpha = (r_PP + r_AP + r_SP + r_PC + r_AC + r_SC) / R0
		Gamma = (r_PC + r_AC + r_SC) / (r_PP + r_SP + r_AP + r_PC + r_AC + r_SC)
		# cat(v," ",R0," ",Theta," ",Omega," ",Alpha,"\n")
		rtnval = 	((Theta - target_Theta)/target_Theta)^2 + 
				((R0 - target_R0)/target_R0)^2 +
				((Omega - target_Omega)/target_Omega)^2  + 
				((Alpha - target_Alpha)/target_Alpha)^2  +
				((Gamma - target_Gamma)/target_Gamma)^2
	} 
	rtnval 
}

no_params = 5
output = array(0,c(no_samples,1+2*no_params))
output = as.data.frame(output)
names(output) <- c("Index","R0","Theta","Omega","Alpha","Gamma","beta","h_pre_post","h_asy_pre","h_pg_hh","h_comm_pg")
output$Index = 1:no_samples

for (i in 1:no_params) {
	if (no_samples==1) output[,1+i] = single_set[i]
	else {
		tmp_H_basev = hyper_vector(no_samples,hc_lb[i],hc_ub[i],hc_log[i])
		output[,1+i] = tmp_H_basev
	}
}

rtnvals = array(0,c(no_params))
current = array(0,c(length(upper_bounds)))

for (j in 1:no_samples) {

	target_R0 = output$R0[j]
	target_Theta = output$Theta[j]
	target_Omega = output$Omega[j]
	target_Alpha = output$Alpha[j]
	target_Gamma = output$Gamma[j]

	min_div = 1000

	for (i in 1:no_iterations) {
		for (k in 1:length(lower_bounds)) current[k]=runif(1,min=lower_bounds[k],max=upper_bounds[k])[1]	
		result = optim(current,reparam,method="L-BFGS-B",lower=opt_lower_bounds,control=list(parscale=scale_factors))
		## result = optim(current,reparam,method="SANN")
		val = result$value[1]
		if (val < min_div) {
			min_div = val
			rtnvals = result$par
		}
		if (val < 10 || i==1) cat(j,"\t",i,"\t",val,"\t",result$par,"\t",rtnvals,"\n")

		flush.console()
	}
	output[j,(1+no_params+1):(1+2*no_params)] = rtnvals
}

output
write.table(output,file=outfile,sep="\t",row.names=FALSE)