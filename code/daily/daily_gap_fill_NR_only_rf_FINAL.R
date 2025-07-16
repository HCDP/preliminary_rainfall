# This code is designed to fill all gaps in the daily rainfall time series.
# Drawing on input files that were previously generated in a companion code.
# Fills all rf stations that had daily RF records with last 2 years

rm(list = ls())
library(dplyr)
library(matrixStats)
library(geosphere)

#set dirs
mainDir <- "/home/hawaii_climate_products_container/preliminary"
codeDir<-paste0(mainDir,"/rainfall/code/source")

#in dirs
inDir <- paste0(mainDir,"/rainfall/data_outputs/tables/station_data/daily/raw_qc/statewide") #Path to the QAQC daily data
activeStaDir <- paste0(mainDir,"/rainfall/data_outputs/tables/rf_station_tracking/lastObs") #last RF obs date per station file dir
gapFilesDir <- paste0(mainDir,"/rainfall/dependencies/daily/gapFilling") #gap relationship files dir

#out dirs
rfFilledDir<-paste0(mainDir,"/rainfall/data_outputs/tables/station_data/daily/partial_filled/statewide") #partial gap filled rf data dir
gapLogDir <- paste0(mainDir,"/rainfall/data_outputs/tables/rf_station_tracking/gapFilling") #gap log dir

#get metadata
meta_url <- "https://raw.githubusercontent.com/ikewai/hawaii_wx_station_mgmt_container/main/Hawaii_Master_Station_Meta.csv"
geog_meta<-read.csv(meta_url, colClasses=c("NESDIS.id"="character"))

#define date
source(paste0(codeDir,"/dataDateFunc.R"))
dataDate<-dataDateMkr() #function for importing/defining date as input or as yesterday
doi<-dataDate #dataDate as date of interest (doi)
file_date<-format(doi,"%Y_%m")

#Read in monthly table
setwd(inDir)
Daily_RF <- read.csv(file = paste0("Statewide_QAQC_Raw_Daily_RF_mm_",file_date,".csv"))
#head(Daily_RF)

#Get the stations that have reported in last 24 months
setwd(activeStaDir)
lastObs<-read.csv("lastRFdayObs.csv")
lastObs$lastRFDate<-as.Date(lastObs$lastRFDate)
lastObs2yr<-lastObs[lastObs$lastRFDate>(dataDate-365.25*2),]
N_Sta <- nrow(lastObs2yr)

#get date cols
dateCols<-grep("X",names(Daily_RF))
minHistCol<-min(dateCols)
maxHistCol<-max(dateCols)

#add stations that reported in last 2 years not already found in table
addSKNs<-lastObs2yr[!lastObs2yr$SKN %in% Daily_RF$SKN,"SKN"]
addDF<-data.frame(matrix(nrow = length(addSKNs), ncol = length(dateCols)+1))
names(addDF)<-names(Daily_RF)[c(1,dateCols)]
addDF$SKN<-addSKNs
addDF<-merge(addDF,geog_meta,by="SKN")
Daily_RF<-rbind(Daily_RF,addDF)
head(Daily_RF)

# Set up Gapfilled output   
RF_Filled <- Daily_RF #duplicate the Input rf file 

#Day of Month Loop
#Loop through each day and fill any missing values
for (D in minHistCol:maxHistCol)  {    #Start minHistCol after metadata ends

  # keep track of day
  print(paste("filling:",gsub("X","",names(Daily_RF)[D])))

  # Cbind SKN and daily RF
    RF_DAY <- cbind(Daily_RF[1],Daily_RF[D])

    # Station Loop 
      for (S in 1:N_Sta) {
        #print(S)
      
    # Get SKN for the TARGET STAION
      TARG_SKN   <- Daily_RF[S,1]  
 
    # Get Value to be filled     
      TARG_Value <- Daily_RF[S,D]
    
    # If Target Value is missing (NA) than proceed
      if(is.na(TARG_Value)) { 
    
    # Call station Input
    # All stations should have an input file with at least distance to other stations.
    # This is a logcial test to make sure there is an Input file based on SKN    
      StaInputTest <- length(list.files(path = gapFilesDir, patt = paste0("SKN_",TARG_SKN,"_Input.csv" ) )) > 0
  
    # If Input file is available 
      if(StaInputTest == TRUE) {   

    # Read in input file 
      StaInfo <- read.csv(paste0(gapFilesDir,"/SKN_",TARG_SKN,"_Input.csv" ))
    
    # Merge Input data and Rainfall Data 
      INFO_RF <- merge(x = StaInfo,y = RF_DAY,by = "SKN")
       
    # Here we choose between two gap-fill techniques based on the result of this test         
    
    # First check to see if the NR Method can be used.
    # Organize the INPUt by Spearmann Rank correlation 
      NR_Input1 <- INFO_RF[order(-INFO_RF$Spear),] #Organize by strongest correlation
    # Sort By NA for rainfall  
      NR_Input2 <- NR_Input1[order(is.na(NR_Input1[minHistCol])),]
    
    # Multipley Ratio by RF Day 
      NR_RAT <- NR_Input2[10] * NR_Input2[minHistCol]
    # Add to DF
      NR_Input2[minHistCol+1]<-NR_RAT
     
    } #End Station Loop for a given day
    
   #Conditions necessary to run NR(most of these alredy adressed)
    # Must have input file
    # Must have a Spearman rank correlation > 0.6 
    # Must have Rainfall   
    # Must have 1100 overlapping points 
    # Ratio must be > 0.3 and less than 2.7
      if(StaInputTest == TRUE && !is.na(NR_Input2[1,minHistCol+1])) {
   
     # Predicted rainfll is the average of predictions at 3 predictor stations
      NR_PRED  <- mean(NR_Input2[1:3,minHistCol+1],na.rm=TRUE) 
      
    # Add precautionaly measure to avoid artificial RF
    # If the highest correlated staion is zero RF and st target station to zero  
      if(NR_Input2[1,minHistCol] < 0.15){NR_RF <- 0}
      
    # Add NR predicted value to Daily time series
      RF_Filled[S,D] <- NR_PRED 
      
    
    } #END NR Fill Statement
  
##############################################
############################################      
    
  # If there is no input file or NR method is not available then do this
  # Note here should always be an input file but in case there isn't then the code won't break down
  # This part of the code will use the distances in the available month 
   else { 
     
     RF_Filled[S,D] <- NA
   
    }  # End else (Choice of Gap-Filling method Fill or no Fill)
   
    }  # End Missing value test (if no missing than move on)
    
    }  # End Station Loop

    }  # End Day of Month Loop 


#remove station rows where all date cols are NA
if(length(dateCols)){}else{}
datesFill<-as.data.frame(RF_Filled[,dateCols]) #only date cols
names(datesFill)<-names(RF_Filled)[dateCols] #rename date cols with date col names
RF_Filled<-RF_Filled[rowSums(is.na(datesFill)) != ncol(datesFill), ] #remove rows where all date cols are NA

#Write (or overwrite) Final Gap-filled dataset
setwd(rfFilledDir) #set rainfall output wd
parFillRFName <- paste0("Statewide_Partial_Filled_Daily_RF_mm_",file_date,".csv")
ifelse(file.exists(parFillRFName),paste(parFillRFName,"file updated!"),paste(parFillRFName,"file written!"))
write.csv(RF_Filled,parFillRFName,row.names = FALSE)

#make and save statistics to say how much of the data was filled. 
#unfilled
naUnfilled<- sum(is.na(Daily_RF[,minHistCol:maxHistCol]))
rfUnfilled<- sum(!is.na(Daily_RF[,minHistCol:maxHistCol]))
#filled
naParFilled<- sum(is.na(RF_Filled[,minHistCol:maxHistCol]))
rfParFilled<- sum(!is.na(RF_Filled[,minHistCol:maxHistCol]))

#Percent NA after fill 
unfilledNAPer<-round(naUnfilled/(rfUnfilled+naUnfilled),3)*100
filledNAPer<-round(naParFilled/(rfParFilled+naParFilled),3)*100
perRFchange<-round((rfParFilled-rfUnfilled)/rfUnfilled,3)*100
perNAchange<-round((naParFilled-naUnfilled)/naUnfilled,3)*100  
rf_Gap_Log<-data.frame(date=Sys.Date(),naUnfilled,rfUnfilled,naParFilled,rfParFilled,unfilledNAPer,filledNAPer,perRFchange,perNAchange)

#write or append file
setwd(gapLogDir)
gap_log_month_filename<-paste0(file_date,"_gap_fill_log_daily_rf.csv")#dynamic file name that includes month year so when month is done new file is written

#conditional statement that adds obs of per day station counts
if(file.exists(gap_log_month_filename)){
  write.table(rf_Gap_Log,gap_log_month_filename, row.names=F,sep = ",",col.names = F, append = T)
  print(paste(gap_log_month_filename,"daily gap log appended!"))
}else{
  write.csv(rf_Gap_Log,gap_log_month_filename, row.names=F)
  print(paste(gap_log_month_filename,"daily gap log written!"))
}

print("gap fill PAU!")
#code end

    