# Set the working directory to the project directory
setwd("../H1N1pdm")

# Ordering is important bacuse of rm(list=ls(all=TRUE)) commands
source("./linelists.R")
source("./fig_linelist.R")
source("./ventilation.R")
source("./fig_baseline.R")
source("./fig_sens.R")
source("./fig_WHO.R")
source("./fig_FHB.R")
source("./fig_early_us_mex.R")
