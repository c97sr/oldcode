
##########################################
#					 #
#  functions for reformatting the data   #
#					 #
##########################################

######################################
# functions for 1 household as 1 entry

#1.3.1.1 function for get 1 household things
getallhousehold_thing <-function(hh, namelist, entryname){
    
    result_data <- hh[entryname]
    names(result_data) <- namelist
    
    hhID <- as.character(hh$q2_h_no) 
    scrID <- as.character(hh$Record_ID)
    result_data <- cbind (scrID, result_data)
    result_data <- cbind (hhID, result_data)

    return (result_data)
}
#end


############################################
#functions for 1 household member as 1 entry

# 1.3.2.2 function for get 1 member things
get1member_thing <-function(hh, household, namelist, entryname){

if (is.na(hh[household,]$q2_1)){ onelength <-9}
else {onelength <- hh[household,]$q2_1} # check the household size by subject self reported on Q2

#one household data
 vectorlen <- length(namelist)+3
 data<- matrix(ncol = vectorlen, nrow=onelength) #set matrix dimension
 
 housename <- c("hhID","scrID","member")
 colnames(data) <- c(housename, namelist) #set matrix names
  

# add ID column
for (i in 1:onelength){
  data[i,1] <- as.character(hh[household,]$q2_h_no) #add household number
     
  data[i,2]  <-  as.character(hh[household,]$Record_ID)  #add screening number    
  data[i,3] <- (i-1)  #add member column
}
  

for (i in 1:onelength){ # i is the number of household member
  member <- i -1

  # local function to input the data
  input_to_matrix <-function( hh, household, name, i, j, data){
	entry_num <-get_entry_num(hh, name)
	data[i, j] <-as.character (hh[household, entry_num])
	return(data)
  }
  #end  local function
   
  for (j in 1: length(entryname)){
  
  data <-  input_to_matrix(hh, household, paste(entryname[j], member, sep="") , i, j+3, data) 
  }
} 
 return(data)

}
#end 


# 1.3.2.1 function for output all member things
getallmember_thing <-function(hh,name,entry){
   output <- get1member_thing(hh,1,name,entry)
   
   for (x in 2:dim(hh)[1]){
   addmatrix <-  get1member_thing(hh,x,name,entry)
   output<- rbind(output, addmatrix)
   }

 output <- as.data.frame(output)
 return(output)
}
#end 



####################################################
#functions for 1 houshold member in 1 day as 1 entry

# 1.3.3.2 function for get 1 member 1 day things
get1day_thing <-function(hh, household, namelist, entryname){
#post:
if (is.na(hh[household,]$q2_1)){ onelength <-9}
else {onelength <- hh[household,]$q2_1} # check the household size by subject self reported on Q2

#one household data
 vectorlen <- length(namelist)+4
 data<- matrix(ncol = vectorlen, nrow=onelength*10) #set matrix dimension
 
 housename <- c("hhID","scrID","member","day")
 colnames(data) <- c(housename, namelist) #set matrix names
  

# add ID column
for (i in 1:(onelength*10)){
 data[i,1]  <-  as.character(hh[household,]$q2_h_no) # add household number
 data[i,2]  <-  as.character(hh[household,]$Record_ID)  #add screening number
}


#add member and day column
count <-1
for (i in 1:onelength){
  for (j in 1:10){
     data[ count ,3 ] <- (i-1)		#add member column
     count <-count + 1
     data[ ( (i-1) *10)+j ,4 ] <- (j-1)    # add day column
   }  
 }

# add entries

for (i in 1:onelength){	# i is the number of household member
 for (whatday in 1:10){ # whatday is the day
   member <- i -1 
   day <- whatday-1

   #local function to input the data entry by entry
  input_to_matrix <-function( hh, household, name, i, day, j, data){
	entry_num <-get_entry_num(hh, name)
	
	data[(i-1)*10 + day, j] <-as.character (hh[household, entry_num])
	return(data)
  }
  #end function


 for (j in 1: length(entryname)){
   thename <- paste("q3_m", member,"d",day,"_", entryname[j], sep="")  # get the database name of the entry 
 
  data <-  input_to_matrix(hh, household, thename , i, whatday, j+4, data) #extract the database's name data out
   }
  }
} 
 return(data)

}
#end



# 1.3.3.1 function for getting all mamber all day things out
getallday_thing <-function(hh,name,entry){
   output <- get1day_thing(hh,1,name,entry)
   
   for (x in 2:dim(hh)[1]){
   addmatrix <-  get1day_thing(hh,x,name,entry)
   output<- rbind(output, addmatrix)
   }

 output <- as.data.frame(output)
 return(output)
}
#end




#functions for all 

#1.3.4.1 extract the data that use comma separated data, return a vector list with boolean entries

#extract1_list_data <- function(datalist, dataname, datavar){
#
#  onelength <- length(dataname)
#  vector_data <-vector(mode= "numeric", length=onelength)
#
#  names(vector_data) <- dataname
#  
# x <- as.character(datalist)
# 
# if (is.na(datalist) | datalist=="")  { vector_data <- vector_data }
#
# else{
#   y <-strsplit(x,split=",")     # split those coma separated symptoms
#   y <-as.numeric(y[[1]])   
#
# for (i in 1:length(y)){
#   for (j in 1:length(datavar)) {
#   if (y[i]  == datavar[j] & !is.na(y[i])) vector_data[j] <-1    
#   }
# }
#}
#return(vector_data)
#}
#end

extract1_list_data <- function(datalist, dataname, datavar){

  onelength <- length(dataname)
  vector_data <- rep(NA,onelength)

  names(vector_data) <- dataname
  
  x <- as.character(datalist)
 
   y <-strsplit(x,split=",")     
   y <-as.numeric(y[[1]])  
   
   if(!is.na(y)){
      vector_data[1:onelength] <-0 
        if(y!=0){
          for (i in 1:length(y)){
          vector_data[y[i]-2] <-1    
	  }
        }
    }
return(vector_data)
}
#end

#1.3.4.2 function for adding the list data in all entry and concatenate at the back 
addall_list_data <-function(dataset, list_entry, list_name, list_var){
  all_sep <- extract1_list_data(list_entry[1], list_name, list_var)
 
  for (x in 2:dim(dataset)[1]) {
  temp_sep <- extract1_list_data(list_entry[x], list_name, list_var)

  all_sep <- rbind(all_sep, temp_sep)
 
  }
data_output <- cbind(dataset, all_sep)
return (data_output)

}
#end


#1.3.4.3 function for deleting missing entries in the database after data extraction
del_miss <- function(data, datatype){


  tocount <- 2
  if (datatype == "household") {tocount<- 1}
  else if (datatype == "member") {tocount<- 4}
  else if (datatype == "day") {tocount<- 5}

  count <- 0
  for (i in 1:dim(data)[1] ){
    for (j in tocount:dim(data)[2] ){
      if (is.na(data[i,j]) | data[i,j]=="")  { count <- count + 1 }
    }
    data$na_count[i] <- count
   
    count <- 0
  } 
 
  cleandata <- data[data$na_count < (dim(data)[2] - tocount),] #remove entries with all NAs
  cleandata <- cleandata[!is.na(cleandata$hhID),]  #remove entries with no hhID
  num <- get_entry_num(cleandata,"na_count" ) #delte the dummy variable for delete the NAs
  cleandata <- cleandata[-c(num)]

  return(cleandata)
 
}  
#end




#1.3.4.4 change dob to age
dob_to_age <-function (dob){

  require(chron)
  today <- "10/10/07"

dob <- dates(as.character(dob), format= "m/d/y")
age <-round((dates(today, format= "m/d/y") - dob)/365,0)
if (age <0) age <-age +100
return(age)
}
#end

#1.3.4.5 change the 2 into 0 in those entries
change2to0 <-function(entry){
  change_entry <- entry
  change_entry <- as.character(change_entry)
  change_entry[!is.na(change_entry) & change_entry == "2" ] <-0
  return(change_entry)
 }

#1.3.4.6 change NA to 0 in those entries
changeNAto0 <-function(entry){
  change_entry <- entry
  change_entry <- as.character(change_entry)
  change_entry[is.na(change_entry) ] <-0
  return(change_entry)
 }





#other functions not used but maybe useful in future

# 1. old delete the missing (NA) data  in the dataset for 1 household member type
del_missing <- function(data){


  count <- 0
  for (i in 1:dim(data)[1] ){
    for (j in 4:dim(data)[2] ){
      if (is.na(data[i,j]))  { count <- count +1 }
    }
    data$na_count[i] <- count
   
    count <- 0
  } 
 
  cleandata <- data[data$na_count < (dim(data)[2] - 4),] #remove entries with all NAs
  cleandata <- cleandata[!is.na(cleandata$hhID),]  #remove entries with no hhID

  return(cleandata)
 
}  
#end





#2. function for getting the smaller dataset from which entry to which entry
get_smaller_data <-function( data, begin_name, end_name){
  
  #function for checking the row number of an entry
  get_entry_num<-function(data ,entry_name){
	return (match(entry_name, names(data)))
  }

  begin_num <- get_entry_num(data, begin_name)
  end_num <- get_entry_num(data, end_name)
  small_data <- data[,begin_num:end_num]
  data_with_ID <- insert_ID(data, small_data)

  return (data_with_ID)
}
#end
# 3. function for getting the entry number by inputting the entry name
 get_entry_num<-function(data ,entry_name){
	return (match(entry_name, names(data)))
  }
#end
#4. insert the household and screening number into the smaller datasets
insert_ID <- function(big_data, small_data) {
	ID <- big_data["ID"]
	h_no <- big_data["q2_h_no"]
	temp <- cbind(h_no,small_data)
	temp <- cbind(ID, temp)
	return (temp)
}
#end

