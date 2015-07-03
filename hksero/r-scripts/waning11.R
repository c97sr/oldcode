# R file to summerize lab results
# Assumes that thew R object in tmp_pp_data.Rdata are up to date
# Make sure that "severity.R" has been run successfully if you are unsure
# If needed: setwd("~/Dropbox/svneclipse/hksero/r-scripts")

rm(list=ls(all=TRUE))
source("../../H1N1pdm/funcs.R")
source("funcs_hksero.R")

# Load up the object pp_data
load("../anon_data/tmp_pp_data.Rdata")
studyTab <- pp_data$base

# Load up the most recent lab results
labTab <- read.csv("../anon_data/1_yr_after_titration_CA4_mahen_sr.csv",stringsAsFactors=FALSE)

# Make a new table for data analysis
labTab$pre_pandemic <- lapply(labTab$pre_pandemic,convertRawTitre) 
labTab$post_pandemic <- lapply(labTab$post_pandemic,convertRawTitre) 
labTab$one_yr_titer <- lapply(labTab$one_yr_titer,convertRawTitre)

labTab$age <- studyTab$age[studyTab$ind_id %in% labTab$Original_id]
(studyTab[studyTab$ind_id %in% labTab$Original_id,])$age

pdf("../outputs/delta_neuts.pdf")
plot(labTab$age,as.numeric(labTab$one_yr_titer) - as.numeric(labTab$post_pandemic))
dev.off()

# Playing around with database extraction functions
x <- mainTab[mainTab$ind_id %in% c("S090586-3","S090586-1"),]

