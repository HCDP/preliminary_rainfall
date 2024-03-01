#this  code grabs 24 hr rainfall total from NWS site from the midnight data upload

rm(list = ls())#remove all objects in R

#packages 
require(plyr)

#settings
options(warn=-1)#suppress warnings for session
Sys.setenv(TZ='Pacific/Honolulu') #set TZ to honolulu
print(paste("nws rr5 daily rf agg:",Sys.time()))#for cron log

#set MAIN DIR
mainDir <- "/home/hawaii_climate_products_container/preliminary"
codeDir<-paste0(mainDir,"/rainfall/code/source")

#define dates
source(paste0(codeDir,"/dataDateFunc.R"))
dataDate<-dataDateMkr() #function for importing/defining date as input or as yesterday
currentDate<-dataDate #dataDate as currentDate

#dirs
parse_hrly_wd<-paste0(mainDir,"/data_aqs/data_outputs/nws_rr5/parse")
agg_daily_wd<-paste0(mainDir,"/rainfall/working_data/nws_rr5")

#read NWS parsed table from dev server OLD
# setwd(parse_hrly_wd)
# nws_filename<-paste0(format((currentDate),"%Y%m%d"),"_nwsrr5_parsed.csv") #dynamic file name that includes date
# nws_hrly_data_final<-read.csv(nws_filename)

#read NWS parsed table from ikewai data portal
ikeUrl<-"https://ikeauth.its.hawaii.edu/files/v2/download/public/system/ikewai-annotated-data/HCDP/workflow_data/preliminary_test" #url
nws_filename<-paste0(format((currentDate),"%Y%m%d"),"_nwsrr5_parsed.csv") #dynamic file name that includes date
nws_hrly_data_final<-read.csv(paste0(ikeUrl,"/data_aqs/data_outputs/nws_rr5/parse/",nws_filename))
#head(nws_hrly_data_final)

#agg to daily and removal partial daily obs
nws_hrly_data_final$date<-as.Date(nws_hrly_data_final$date)#make date a date
nws_hrly_data_final_cc<-nws_hrly_data_final[complete.cases(nws_hrly_data_final),]#remove NA obs
nws_24hr_agg<-ddply(nws_hrly_data_final_cc, .(nwsli,NWS_name, date), summarize, prec_mm_24hr = sum(prec_mm_1hr), hour_count=length(nwsli)) #sum by station & day and count hourly obs per station & day
nws_24hr_agg_final<-nws_24hr_agg[nws_24hr_agg$hour_count >= 23,] #subset by stations with at least 23 hourly obs (ie: only 1 missing data)
nws_24hr_agg_final$nwsli<-as.character(nws_24hr_agg_final$nwsli) #cast nws id to character
airports<-c("HJR","HKO","HLI","HMK","HNL","HNY","HOG","HTO")
nws_24hr_agg_final[nws_24hr_agg_final$nwsli %in% airports,"nwsli"]<-paste0("P",nws_24hr_agg_final[nws_24hr_agg_final$nwsli %in% airports,"nwsli"]) #add P to airport IDs to match meta
#head(nws_24hr_agg_final)
#tail(nws_24hr_agg_final)

#write/append to file 
setwd(agg_daily_wd)#set wd to save/append final day aggs
files<-list.files()
rf_month_daily_filename<-paste0(format((currentDate),"%Y_%m"),"_nws_daily_rf.csv")#dynamic filename that includes month year so when month is done new file is written
if(file.exists(rf_month_daily_filename)){
	write.table(nws_24hr_agg_final,rf_month_daily_filename, row.names=F, sep = ",", col.names = F, append = T)
	print(paste(rf_month_daily_filename,"appended"))
	}else{
	write.csv(nws_24hr_agg_final,rf_month_daily_filename,row.names=F)
	print(paste(rf_month_daily_filename,"written"))
	}
print("PAU!")

# CODE PAU!!!!!
