#
# Figure 1: flowchart
#

# Totally 946 index subjects involved (2 excluded b/c hard copy record missing)
# 198 randomized (meet inclusion criteria), 70 refused to participate
# 128 allocated to control arm, 35 allocated to mask arm, 35 allocated to hand hygiene arm

qv <- read.csv("P:\\GROUP\\NPIpilot\\data\\2007_11_28_clinicdat_h.csv", header=TRUE)
hchar <- read.csv("P:\\GROUP\\NPIpilot\\data\\2007_11_27_housechar_h.csv", header=TRUE)
hc <- read.csv("P:\\GROUP\\NPIpilot\\data\\2007_11_16_home_culture.csv", header=TRUE)
baseflu <- read.csv("P:\\GROUP\\NPIpilot\\data\\2007_11_20_baseflu_m.csv", header=TRUE)

# Step 1: Enrolment

# QuickVue +ve among 198 index subjects
dim(qv[qv$QVres==1 | qv$QVres==2,])[1]          # QV: +ve
dim(qv[qv$QVres==1,])[1]                        # QV: A
dim(qv[ qv$QVres==2,])[1]                       # QV: B

# End of step 1.#####################################################################################################################

# Step 2: Allocation

# Received allocated intervention
control <- hchar[hchar$intervention==1,]              # Control arm
dim(control)[1]
median(control$q1_familysize)
range(control$q1_familysize)
sum(control$q1_familysize) - dim(control)[1]

mask <- hchar[hchar$intervention==2,]                 # Mask arm
dim(mask)[1]
median(mask$q1_familysize)
range(mask$q1_familysize)
sum(mask$q1_familysize) - dim(mask)[1]

hand <- hchar[hchar$intervention==3,]                 # Hand hygiene arm
dim(hand)[1]
median(hand$q1_familysize)
range(hand$q1_familysize)
sum(hand$q1_familysize) - dim(hand)[1]

# End of step 2.#####################################################################################################################

# Step 3: Analysis

hculture <- hc[hc$member==0, c(1:3,6,8:11)]           # Extract home culture results for index subhects

hc0 <- hculture[hculture$visit==0,]                   # Extract home culture results for index subhects (visit 0)
for (j in 1:nrow(hc0)){
     if( (!is.na(hc0$culture[j])&as.character(hc0$culture[j])>0) | (!is.na(hc0$AqPCR[j])&as.character(hc0$AqPCR[j])>0) |
             (!is.na(hc0$AqPCR2[j])&as.character(hc0$AqPCR2[j])>0) | (!is.na(hc0$AqPCR3[j])&as.character(hc0$AqPCR3[j])>0) |
	     (!is.na(hc0$BqPCR[j])&as.character(hc0$BqPCR[j])>0) )
	     hc0$cultpos0[j] <- 1
	     else  hc0$cultpos0[j] <- 0
}

hc1 <- hculture[hculture$visit==1,]                   # Extract home culture results for index subhects (visit 1)
for (j in 1:nrow(hc1)){
     if( (!is.na(hc1$culture[j])&as.character(hc1$culture[j])>0) | (!is.na(hc1$AqPCR[j])&as.character(hc1$AqPCR[j])>0) |
             (!is.na(hc1$AqPCR2[j])&as.character(hc1$AqPCR2[j])>0) | (!is.na(hc1$AqPCR3[j])&as.character(hc1$AqPCR3[j])>0) |
	     (!is.na(hc1$BqPCR[j])&as.character(hc1$BqPCR[j])>0) )
	     hc1$cultpos1[j] <- 1
	     else  hc1$cultpos1[j] <- 0
}

hc01 <- merge(hc0[,c(1,2,9)],hc1[,c(1,9)],by="hhID",all.x=T)               # Define 'baseline'(=1) when V0 & V1 culture -ve
hc01$baseline <- 1*(!is.na(hc01$cultpos0) & hc01$cultpos0==0 & !is.na(hc01$cultpos1) & hc01$cultpos1==0)  
hc_comb01 <- hc01[,c(1,5)]
analysis <- merge(hchar,hc_comb01,by="hhID",all.x=T)


# Analysed clusters:
control_als <- analysis[analysis$baseline==0 & analysis$intervention==1,]
dim(control_als)[1]
median(control_als$q1_familysize)
range(control_als$q1_familysize)
sum(control_als$q1_familysize) - dim(control_als)[1]           # Participants

mask_als <- analysis[analysis$baseline==0 & analysis$intervention==2,]
dim(mask_als)[1]
median(mask_als$q1_familysize)
range(mask_als$q1_familysize)
sum(mask_als$q1_familysize) - dim(mask_als)[1]                 # Participants

hand_als <- analysis[analysis$baseline==0 & analysis$intervention==3,]
dim(hand_als)[1]
median(hand_als$q1_familysize)
range(hand_als$q1_familysize)
sum(hand_als$q1_familysize) - dim(hand_als)[1]                 # Participants


# End of step 3.#####################################################################################################################
