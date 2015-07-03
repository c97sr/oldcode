# Create a new .csv source file that indicates whether a specific household was recruited through random dialing or from Swineflu study
rm(list=ls(all=TRUE))
options(error=NULL)

# Read input files
lab_results2_raw <- read.csv("M:\\raw_data\\Neutralization_Test_Result_of_Steven_sera_(Titration)_May_12_2010.csv")

# Initalize variables
norowlab2 <- nrow(lab_results2_raw)
datafile <- matrix(0, nrow=norowlab2, ncol=8)

# New variables
Base_titre <- array(-1,c(1,norowlab2))
Post_titre <- array(-1,c(1,norowlab2))
test2_ratio <- array(-1,c(1,norowlab2))
Kpos <- array(-1,c(1,norowlab2))

# Main program
for (i in 1:norowlab2) {

        # Extract baseline titre levels
        if ("Negative (<1:10)" %in% lab_results2_raw$Baseline[i]) Base_titre[i] <- 5 else
        if ("1:10" %in% lab_results2_raw$Baseline[i]) Base_titre[i] <- 10 else
        if ("1:20" %in% lab_results2_raw$Baseline[i]) Base_titre[i] <- 20 else
        if ("1:40" %in% lab_results2_raw$Baseline[i]) Base_titre[i] <- 40 else
        if ("1:80" %in% lab_results2_raw$Baseline[i]) Base_titre[i] <- 80 else
        if ("1:160" %in% lab_results2_raw$Baseline[i]) Base_titre[i] <- 160 else
        if ("1:320" %in% lab_results2_raw$Baseline[i]) Base_titre[i] <- 320 else
        if ("1:640" %in% lab_results2_raw$Baseline[i]) Base_titre[i] <- 640 else
        print('Something is missing in Baseline titre')
        
        # Extract post-season titre levels
        if ("Negative (<1:10)" %in% lab_results2_raw$Post[i]) Post_titre[i] <- 5 else
        if ("1:10" %in% lab_results2_raw$Post[i]) Post_titre[i] <- 10 else
        if ("1:20" %in% lab_results2_raw$Post[i]) Post_titre[i] <- 20 else
        if ("1:40" %in% lab_results2_raw$Post[i]) Post_titre[i] <- 40 else
        if ("1:80" %in% lab_results2_raw$Post[i]) Post_titre[i] <- 80 else
        if ("1:160" %in% lab_results2_raw$Post[i]) Post_titre[i] <- 160 else
        if ("1:320" %in% lab_results2_raw$Post[i]) Post_titre[i] <- 320 else
        if ("1:640" %in% lab_results2_raw$Post[i]) Post_titre[i] <- 640 else
        print('Something is missing in post-season titre')  
}
        
        
for (j in 1:norowlab2) {
    if (Base_titre[j]==5 && Post_titre[j]!=5) test2_ratio[j] <- Post_titre[j]/ 10 else
    if (Base_titre[j]!=5 && Post_titre[j]==5) test2_ratio[j] <- 10/ Base_titre[j] else
    if (Base_titre[j]==5 && Post_titre[j]==5) test2_ratio[j] <- 10/10 else
    test2_ratio[j] <-  Post_titre[j] / Base_titre[j]
        
    if (test2_ratio[j]>=4) Kpos[j] <- 1 else
            Kpos[j] <- 0
}
            
datafile[,1] <- as.character(lab_results2_raw$ind_id)
datafile[,2] <- as.character(lab_results2_raw$ind_id2)
datafile[,3] <- as.character(lab_results2_raw$Baseline)
datafile[,4] <- as.character(lab_results2_raw$Post)
datafile[,5] <- Base_titre
datafile[,6] <- Post_titre
datafile[,7] <- test2_ratio
datafile[,8] <- Kpos
colnames(datafile) <- c("ind_id","ind_id2","Baseline","Post","Base_titre","Post_titre","Ratio","Kpos")

# Save subdatasets
write.csv(datafile, file="M:\\svn\\anon_data\\lab_results2.csv")
