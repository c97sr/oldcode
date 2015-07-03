# This file needs to be adapted to work with the updated version of the amoy proceedures

rm(list=ls(all=TRUE))
# below if needed
# setwd("")

source("time_since_functions.r")

# infile = "D:\\files\\projects\\sars\\amoy\\data\\flu_sim_input.in"
# outfile = "D:\\files\\projects\\sars\\amoy\\data\\flu_output.out"

infile = "flu_sim_input.in"
outfile = "flu_output.out"

amPs <- data.frame(
	q_ptp = 0.05,
	mu_ptp = 4,
	alpha_ptp = 3,
	t_off_ptp = 0,
	q_seed = 0, 				
	mu_seed = 1,
	alpha_seed = 100,
	onset_ave_dur = 4.65,
	T_B = 1000000
)

op <- amoy_sim(infile,amPs,outfile,seed=2948576)

