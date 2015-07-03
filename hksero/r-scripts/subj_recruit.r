# Create a new .csv source file that indicates whether a specific household was recruited through random dialing or from Swineflu study
rm(list=ls(all=TRUE))
options(error=NULL)

# Read input files
id_raw <- read.csv("M:\\raw_data\\serosurvey_id_20100210.csv")
telsource_raw <- read.csv("M:\\raw_data\\Blood09_namecontact_POP_2010may26.csv")

# Initalize variables
norowid <- dim(id_raw)[1]
datafile <- matrix(0, nrow=norowid, ncol=13)
telmatch1 <- array(FALSE, c(1,norowid))
telmatch2 <- array(FALSE, c(1,norowid))
telmatch3 <- array(FALSE, c(1,norowid))
telmatch4 <- array(FALSE, c(1,norowid))
telmatch <- array(FALSE, c(1,norowid))
telsource1 <- array(-1, c(1,norowid))
telsource2 <- array(-1, c(1,norowid))
telsource3 <- array(-1, c(1,norowid))
telsource4 <- array(-1, c(1,norowid))
telsource <- array(-1, c(1,norowid))

# Main program
for (i in 1:norowid) {

        # Telephone number matched or not
        if (id_raw$home_phone[i] %in% telsource_raw$tel_no) telmatch1[i] <- TRUE
        else if (id_raw$home_phone[i] %in% telsource_raw$contact) telmatch2[i] <- TRUE
        else if (id_raw$mobile_phone[i] %in% telsource_raw$tel_no) telmatch3[i] <- TRUE
        else if (id_raw$mobile_phone[i] %in% telsource_raw$contact) telmatch4[i] <- TRUE        
     
        if (telmatch1[i]==TRUE || telmatch2[i]==TRUE || telmatch3[i]==TRUE || telmatch4[i]==TRUE) telmatch[i]<-TRUE

        # Source of telephone
        if (telmatch1[i]) {           
            idrow1 <- match(id_raw$home_phone[i],telsource_raw$tel_no)
            if (is.na(idrow1)) print('There should not be a NA here.')
            telsource1[i] <- telsource_raw$tel_source[idrow1]
            telsource[i] <- telsource_raw$tel_source[idrow1]
        }
        
        if (telmatch2[i]) {
            idrow2 <- match(id_raw$home_phone[i],telsource_raw$contact)
            if (is.na(idrow2)) print('There should not be a NA here.')
            telsource2[i] <- telsource_raw$tel_source[idrow2]
            telsource[i] <- telsource_raw$tel_source[idrow2]
        }
        
        if (telmatch3[i]) {
            idrow3 <- match(id_raw$mobile_phone[i],telsource_raw$tel_no)
            if (is.na(idrow3)) print('There should not be a NA here.')
            telsource3[i] <- telsource_raw$tel_source[idrow3]
            telsource[i] <- telsource_raw$tel_source[idrow3]
        }
        
        if (telmatch4[i]) {
            idrow4 <- match(id_raw$mobile_phone[i],telsource_raw$contact)
            if (is.na(idrow4)) print('There should not be a NA here.')
            telsource4[i] <- telsource_raw$tel_source[idrow4]
            telsource[i] <- telsource_raw$tel_source[idrow4]
        }
        
        if (telmatch[i] == FALSE) {
        if (id_raw$ind_id[i] == "S090023-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090037-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090045-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090072-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090077-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090080-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090081-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090085-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090086-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090087-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090089-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090090-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090092-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090093-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090094-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090095-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090096-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090097-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090098-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090099-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090100-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090101-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090102-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090103-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090104-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090105-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090106-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090107-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090109-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090110-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090111-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090112-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090115-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090116-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090117-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090118-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090119-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090120-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090121-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090122-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090124-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090125-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090126-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090127-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090128-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090129-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090130-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090131-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090132-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090133-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090134-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090135-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090137-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090138-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090139-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090140-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090141-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090142-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090143-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090144-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090145-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090146-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090159-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090178-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090182-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090184-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090185-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090190-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090191-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090192-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090193-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090197-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090198-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090200-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090202-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090203-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090205-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090206-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090207-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090209-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090212-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090213-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090214-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090215-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090219-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090221-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090224-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090226-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090228-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090229-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090231-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090232-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090233-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090234-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090237-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090238-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090239-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090244-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090249-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090254-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090255-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090259-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090266-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090267-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090272-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090273-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090274-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090284-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090286-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090288-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090289-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090290-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090291-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090294-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090295-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090296-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090298-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090299-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090318-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090324-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090369-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090411-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090412-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090510-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090518-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090553-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090554-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090555-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090556-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090558-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090559-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090560-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090564-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090566-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090567-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090569-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090570-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090571-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090572-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090573-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090575-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090576-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090577-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090579-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090581-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090582-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090584-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} else
        if (id_raw$ind_id[i] == "S090589-0") {telmatch1[i] <- TRUE; telmatch[i] <- TRUE} 
        }
        
        if (telsource[i] == -1) {
        if (id_raw$ind_id[i] == "S090023-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090037-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090045-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090072-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090077-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090080-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090081-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090085-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090086-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090087-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090089-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090090-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090092-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090093-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090094-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090095-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090096-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090097-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090098-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090099-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090100-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090101-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090102-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090103-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090104-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090105-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090106-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090107-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090109-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090110-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090111-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090112-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090115-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090116-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090117-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090118-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090119-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090120-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090121-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090122-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090124-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090125-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090126-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090127-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090128-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090129-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090130-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090131-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090132-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090133-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090134-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090135-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090137-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090138-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090139-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090140-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090141-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090142-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090143-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090144-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090145-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090146-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090159-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090178-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090182-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090184-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090185-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090190-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090191-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090192-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090193-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090197-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090198-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090200-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090202-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090203-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090205-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090206-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090207-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090209-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090212-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090213-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090214-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090215-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090219-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090221-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090224-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090226-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090228-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090229-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090231-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090232-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090233-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090234-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090237-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090238-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090239-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090244-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090249-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090254-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090255-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090259-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090266-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090267-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090272-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090273-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090274-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090278-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090284-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090286-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090288-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090289-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090290-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090291-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090294-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090295-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090296-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090298-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090299-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090318-0") {telsource1[i] <- 2; telsource[i] <- 2} else
        if (id_raw$ind_id[i] == "S090324-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090369-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090411-0") {telsource1[i] <- 2; telsource[i] <- 2} else
        if (id_raw$ind_id[i] == "S090412-0") {telsource1[i] <- 2; telsource[i] <- 2} else
        if (id_raw$ind_id[i] == "S090510-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090518-0") {telsource1[i] <- 2; telsource[i] <- 2} else
        if (id_raw$ind_id[i] == "S090553-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090554-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090555-0") {telsource1[i] <- 2; telsource[i] <- 2} else
        if (id_raw$ind_id[i] == "S090556-0") {telsource1[i] <- 2; telsource[i] <- 2} else
        if (id_raw$ind_id[i] == "S090558-0") {telsource1[i] <- 2; telsource[i] <- 2} else
        if (id_raw$ind_id[i] == "S090559-0") {telsource1[i] <- 2; telsource[i] <- 2} else
        if (id_raw$ind_id[i] == "S090560-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090564-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090566-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090567-0") {telsource1[i] <- 2; telsource[i] <- 2} else
        if (id_raw$ind_id[i] == "S090569-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090589-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090570-0") {telsource1[i] <- 2; telsource[i] <- 2} else
        if (id_raw$ind_id[i] == "S090571-0") {telsource1[i] <- 2; telsource[i] <- 2} else
        if (id_raw$ind_id[i] == "S090572-0") {telsource1[i] <- 2; telsource[i] <- 2} else
        if (id_raw$ind_id[i] == "S090573-0") {telsource1[i] <- 2; telsource[i] <- 2} else
        if (id_raw$ind_id[i] == "S090575-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090576-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090577-0") {telsource1[i] <- 2; telsource[i] <- 2} else
        if (id_raw$ind_id[i] == "S090579-0") {telsource1[i] <- 1; telsource[i] <- 1} else
        if (id_raw$ind_id[i] == "S090581-0") {telsource1[i] <- 2; telsource[i] <- 2} else
        if (id_raw$ind_id[i] == "S090582-0") {telsource1[i] <- 2; telsource[i] <- 2} else
        if (id_raw$ind_id[i] == "S090584-0") {telsource1[i] <- 2; telsource[i] <- 2}
        }
}

datafile[,1] <- as.character(id_raw$ind_id)
datafile[,2] <- id_raw$hh_id
datafile[,3] <- id_raw$hh_index
datafile[,4] <- telmatch1
datafile[,5] <- telmatch2
datafile[,6] <- telmatch3
datafile[,7] <- telmatch4
datafile[,8] <- telmatch
datafile[,9] <- telsource1
datafile[,10] <- telsource2
datafile[,11] <- telsource3
datafile[,12] <- telsource4
datafile[,13] <- telsource
colnames(datafile) <- c("ind_id","hh_id","hh_index","Match1","Match2","Match3","Match4","Match","Source1","Source2","Source3","Source4","Source")

# Save subdatasets
write.csv(datafile, file="M:\\svn\\anon_data\\subj_recruit.csv")
