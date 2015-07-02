# A script to generate the required input files for
# the small scale metapopulation model of transmission in 
# southern China

# Extracts population from landscan
# Generates a list of population age groups for each population
# Generates a table of connectivities between populations

# Minimum requirements for this file are that it
# - can generate a population density file based on landscan and max / min values

# Set the local directory if required setwd("/Users/sriley/Dropbox/svneclipse/globalflu/r-scripts")

# Clear all pbjects and set the error behavior
rm(list=ls(all=TRUE))
options(error=NULL)

require("sp")
# Load up libraries
source("../../R Source/SR_Misc.r")

# Set task toggles
# The section that needs rgdal has only been run on an XP VM
remakepopgrid=TRUE
# Declare some overall parameters
strFnRoot <- "../working/gz_urb_rur"

# Load up the popgrid
fnPopgridData <- "./popgridtmp.Rdata"
if (remakepopgrid) {
	require("rgdal")
	landscanfile <- "Z:Volumes/NO NAME/data/gis/landscan/asia03/asia03/w001001.adf"
	asciifile <- "../working/SmallGZLandScanAscii.txt"
	# These numbers read directly off google earth at the moment 
	ds = 23.1365
	de = 113.2533
	dn = 23.2223
	dw = 113.3520 
	ExcerptLStoASCII(landscanfile,asciifile,de,dw,ds,dn,lsspd=120,lsmn=85,lsmw=34)
	popgrid <- read.asciigrid(asciifile)
	save(popgrid,file=fnPopgridData)
} else {
	load(fnPopgridData)
}

image(popgrid)

# Set up age groups for each population in the grid
# Save as only a csv that will read in as a GSL matrix
# Parameter vector for upper limit of age groups is known 

# Set up a connectivity kernel for each agegroup in the grid
