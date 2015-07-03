# Clear all current objects
rm(list=ls(all=TRUE))

# Function to test if a contact is obese
checkContact <- function(index,aTies,vNoTies,dfChars) {
	rtn <- FALSE
	if (vNoTies[index] > 0)
		for (i in 1:vNoTies[index]) 
			if (dfChars$BMI[aTies[index,i]] > 30) rtn <- TRUE
	rtn
}

# Initialise variables
no_people <- 1000
max_ties <- 20
ave_ties <- 7.5 
no_iterations <- 1
max_time <- 1000
time_step <- 1
trigger_y <- 0.23
final_y <- 0.31
interventions <- TRUE
set.seed(1234)

# Define biological parameters
obPars <- c(
	h_b = 0.004,	# base hazard per year of becoming obese
	h_r = 1.5 		# relative hazard per obese contact
)

ties <- array(NA,c(no_people,max_ties))
ties_each_person <- array(0,c(no_people))
bmi <- array(FALSE,c(no_people))

# Set up a random bidirectional network with no_people nodes 
# of average degree ave_ties. 
# Degree distribution is truncated poisson.
required_ties <- round(no_people * ave_ties)
while (required_ties > 0) {
	first <- ceiling(runif(1)*no_people)
	second <- ceiling(runif(1)*no_people)
	if (!(	is.element(c(first),ties[second,])	 	||
			is.element(c(second),ties[first,]) 		||
			ties_each_person[first] == max_ties 	||
			ties_each_person[second] == max_ties 	||
			first == second		)) {
		ties[first,ties_each_person[first]+1] <- second
		ties[second,ties_each_person[second]+1] <- first
		ties_each_person[first] <- ties_each_person[first]+1
		ties_each_person[second] <- ties_each_person[second]+1
		required_ties <- required_ties - 1
	}
}

# Set up character variables with all BMIs = 22.5
chars <- data.frame("BMI"=array(22.5,c(no_people)))
obese <- vector()
notobese <- 1:no_people

# Initialise simulation loop
for (i in 1:no_iterations) {
	# Seed the epidemic with one individual with BMI=4
	t <- 0
	notstop <- TRUE
	trigger <- FALSE
	while (t <= max_time && notstop == TRUE) {
		newobese <- vector()
		for (j in notobese) {
			if (	checkContact(j,ties,ties_each_person,chars) && 
					!(trigger && interventions)) h_j <- obPars["h_b"]*obPars["h_r"]
			else h_j <- obPars["h_b"]
			# cat(j," ",h_j,"\n")
			if (runif(1)<1-exp(-h_j*time_step)) {
				newobese <- c(newobese,j)
			}
		}
		for (j in newobese) {
			obese <- c(obese,j)
			notobese <- notobese[-match(j,notobese)]
			chars$BMI[j] <- 40
		}		
		y <- length(obese)/no_people
		dy <- length(newobese)/no_people
		t <- t + time_step
		if (y > trigger_y) trigger <- TRUE
		if (y > final_y) notstop <- FALSE
		cat("R\t",i,"\tt\t",t,"\ty\t",y,"\tdy\t",dy,"\n")
		flush.console()	
	}
}