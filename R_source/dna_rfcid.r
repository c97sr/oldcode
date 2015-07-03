# To read in and compare some flu sequences
rm(list=ls(all=TRUE))
require("ape")
# change the directory if required setwd("/Users/sriley/Dropbox/svneclipse/R Source")

# A function to return an auxilliary object for the analysis
setUpDetails <- function(dnaseq) {
	
	# Define some housekeeping variables
	intSeqLength <- length(dnaseq[[1]])
	intNoSequences <- length(dnaseq)
	allDist <- round(intSeqLength*dist.dna(dnaseq,model="raw",pairwise.deletion=FALSE,as.matrix=TRUE))
	# hist(allDist)
	
	# Set up a dataframe for all the sequences
	# For each we should know: Name | Current Study (0|1) | Household | Member | Visit
	dfSeqChars <- as.matrix(cbind(	
					current_study=rep(0,intNoSequences),
					managua_study=rep(0,intNoSequences),
					household=rep(0,intNoSequences),
					member=rep(0,intNoSequences),
					visit=rep(0,intNoSequences))
			)
	
	for (i in 1:intNoSequences) {
		
		# Break the string
		strBits <- (strsplit(names(dnaseq[i]),"/"))[[1]]
		
		# Check for membership of studies
		if (length(strBits) > 1 && strBits[2] == "Hong_Kong") dfSeqChars[i,"current_study"] <- 1
		if (length(strBits) > 1 && strBits[2] == "Managua") dfSeqChars[i,"managua_study"] <- 1
		
		# Assign additional details for current study
		if (dfSeqChars[i,"current_study"] > 0) {
			
			# Split the
			if (length(strBits) < 3) stop("Inconsistancy in setUpDetails: 92837429")
			details <- (strsplit(strBits[3],"m|v"))[[1]]
			if (length(details) != 3) stop("Inconsistancy in setUpDetails: 83737429")
			dfSeqChars[i,"household"] <- as.numeric(details[1])
			dfSeqChars[i,"member"] <- as.numeric(details[2])
			dfSeqChars[i,"visit"] <- as.numeric(details[3])
			
		}
		
	}
	
	list(lengthdna=intSeqLength,noseq=intNoSequences,distances=allDist,chars=dfSeqChars)
	
}

# Set some housekeeping variables
remakeSeqObj = FALSE

# A file of all the full sequences for human H3N2 for which time is available
fn2007AllFullTimed <- "/Users/sriley/Dropbox/projects/influenza/dna_household/anon_data/H3N2_all_concat_154.nex"

if (remakeSeqObj) {

	# Read in the raw data file (with no spaces in the sequence names)
	rawSeqData <- read.nexus.data(fn2007AllFullTimed)
	binSeqData <- as.DNAbin(rawSeqData)
	save(binSeqData,file="binSeqData.Rdata")	
		
	# Alternative example of code reading in fasta files
	# fnHA <- "/Users/sriley/Dropbox/grants/rfcid_dna/sequence_data/ny_2003_clade_A_HA_aligned.fa"
	# HA <- read.dna(fnHA,format="fasta")

} else {

	# Load up the presaved R object
	load("binSeqData.Rdata")
	
}

# Read off the length of the sequences and calculate the distances between each sequence
x <- setUpDetails(binSeqData)

# Construct the frequency tables
max_dist <- max(x$distances)
tabMain <- 	as.matrix(cbind(	host=rep(0,max_dist+1),
								house=rep(0,max_dist+1),
								study=rep(0,max_dist+1),
								study_to_out=rep(0,max_dist+1),
								outgroup=rep(0,max_dist+1)
			))

debugcount <- 0
	
# Cycle though all distinct pair combinations and assign distances
for (i in 1:(x$noseq-1)) {
	
	# Second member of the pair
	for (j in (i+1):x$noseq) {
		
		debugcount <- debugcount + 1
		
		# Distance between pairs
		intCurrentDistance <- x$distances[i,j] + 1
		
		# Test that i and j are members of the current study
		if (	x$chars[i,"current_study"] == 1 && x$chars[j,"current_study"] == 1) {
			
			# Test for the same household
			if (x$chars[i,"household"] == x$chars[j,"household"]) {
				
				# Test for the same individual
				if (x$chars[i,"member"] == x$chars[j,"member"]) {
						
					# Increment distance measure for within host
					tabMain[intCurrentDistance,"host"] <- tabMain[intCurrentDistance,"host"] + 1 
					
				# Close if for same individual open else	
				} else {
					
					# Increment distance measure for within household
					tabMain[intCurrentDistance,"house"] <- tabMain[intCurrentDistance,"house"] + 1 
					
				# Close else for same individual	
				}
					
			
			# Close the if for same household	
			} else {

				# Increment distance measure for within study
				tabMain[intCurrentDistance,"study"] <- tabMain[intCurrentDistance,"study"] + 1 
				
			# Close else for same household	
			}
		
		# Close the if for both current study and open else	
		} else {
			
			# Test if either were in the current study
			if (x$chars[i,"current_study"] == 1 || x$chars[j,"current_study"] == 1) {
				
				# Incremenet the between study count
				tabMain[intCurrentDistance,"study_to_out"] <- tabMain[intCurrentDistance,"study_to_out"] + 1
				
			# Endif of either being in the current study
			} else {
				
				# This case must be that neither are in the current study
				tabMain[intCurrentDistance,"outgroup"] <- tabMain[intCurrentDistance,"outgroup"] + 1
			
			# End else that either are in the current study
			}
			
		}
		
	# Close the j loop
	}

# Close the for i loop
}

if ((x$noseq * x$noseq - x$noseq) / 2 != sum(tabMain)) stop("Problem with the number of pairs")
write.csv(tabMain,file="tabMain.csv")

