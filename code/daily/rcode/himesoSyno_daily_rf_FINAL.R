#this code grabs all raw downloaded mesonet synopic API data for HI and calcs daily RF total from 5 min precip

rm(list = ls())#remove all objects in R

#set global options
options(warn=-1) #suppress warnings for session
Sys.setenv(TZ='Pacific/Honolulu') #set TZ to honolulu
print(paste("mesonet rf daily run:",Sys.time()))#for cron log

#load packages
#install.packages("xts")
require(xts)

#set dirs
mainDir <- "/home/hawaii_climate_products_container/preliminary"
codeDir<-paste0(mainDir,"/rainfall/code/source")
parse_wd<-paste0(mainDir,"/data_aqs/data_outputs/mesonetSynoptic/parse")
agg_daily_wd<-paste0(mainDir,"/rainfall/working_data/hi_mesonet/synoptic")

#define dates
source(paste0(codeDir,"/dataDateFunc.R"))
dataDate<-dataDateMkr() #function for importing/defining date as input or as yesterday
currentDate<-dataDate #dataDate as currentDate

#functions
getmode <- function(v) { #get mode of values
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}#end getmode function
apply.hourly <- function(x, FUN, roundtime = "round", na.rm = TRUE){
  if(!is.xts(x)){
    stop("x must be an xts object")
  }
  
  if(roundtime != "NA"){
    if(roundtime == "round"){
      time(x) <- round.POSIXt(time(x), "hours")
    } else if(roundtime == "trunc"){
      time(x) <- trunc.POSIXt(time(x), "hours")
    } else {
      stop("roundtime must be either round or trunc")
    }
  }
  
  ap <- endpoints(x,'hours')
  if(na.rm){
    period.apply(x,ap,FUN, na.rm = TRUE)
  } else {
    period.apply(x,ap,FUN)
  }
}#end apply.hrly function

# #read meso parsed table from dev server OLD
# setwd(parse_wd)#sever path for parsed meso files
# meso_filename<-paste0(format((currentDate),"%Y%m%d"),"_himeso_synoptic.csv") #dynamic file name that includes date
# all_meso<-read.csv(meso_filename)
# head(all_meso)

#read meso parsed table from ikewai data portal NEW
ikeUrl<-"https://ikeauth.its.hawaii.edu/files/v2/download/public/system/ikewai-annotated-data/HCDP/workflow_data/preliminary_test" #url
meso_filename<-paste0(format((currentDate),"%Y%m%d"),"_himeso_synoptic.csv") #dynamic file name that includes date
all_meso<-read.csv(paste0(ikeUrl,"/data_aqs/data_outputs/mesonetSynoptic/parse/",meso_filename))
#head(all_meso)
#unique(all_meso$var)

#subset precip var, convert inch to mm and convert UTC to HST
all_meso_rf<-subset(all_meso,var=="precip_accum_five_minute_set_1")# subset precip only
all_meso_rf$value<-as.numeric(all_meso_rf$value) #cast value to numeric
all_meso_rf<-all_meso_rf[complete.cases(all_meso_rf),] #remove NA rows
all_meso_rf$Date_Time<-strptime(all_meso_rf$Date_Time, format="%Y-%m-%dT%H:%M:%SZ", tz="UTC")
attr(all_meso_rf$Date_Time,"tzone") <- "Pacific/Honolulu" #convert TZ attribute to HST
all_meso_rf$Date_Time<-(all_meso_rf$Date_Time-36000) #minus 10 hrs for UTC to HST
all_meso_rf$Date_Time<-(all_meso_rf$Date_Time)-1 #minus 1 second to put midnight obs in last day
head(all_meso_rf)
tail(all_meso_rf)

#blank DF to store daily data
meso_daily_rf<-data.frame()

#unique meso stations
stations<-unique(all_meso_rf$Station_ID)

#start daily RF loop
print("daily rf loop started...")
# j<-stations[1]
for(j in stations){
  sta_data<-subset(all_meso_rf,Station_ID==j)
  sta_data_xts<-xts(sta_data$value,order.by=sta_data$Date_Time,unique = TRUE) #make xtended timeseries object
  sta_data_xts<- sta_data_xts[!duplicated(index(sta_data_xts)),] #remove duplicate time obs
  sta_data_hrly_xts<-apply.hourly(sta_data_xts,FUN=sum,roundtime = "trunc")#agg to hourly and truncate hour
  indexTZ(sta_data_hrly_xts) <- "Pacific/Honolulu"
  sta_data_daily_xts<-apply.daily(sta_data_hrly_xts,FUN=sum,na.rm = F) #daily sum of all all lag observations
  obs_ints<-diff(index(sta_data_xts),lag=1) #calculate vector of obs intervals
  obs_int_hr<-getmode(as.numeric(obs_ints, units="hours"))
  obs_int_minutes<-obs_int_hr*60
  obs_per_day<-((1/obs_int_hr)*24)#calculate numbers of obs per day based on obs interval
  sta_per_obs_daily_xts<-as.numeric(apply.daily(sta_data_xts,FUN=length)/obs_per_day)#vec of % percentage of obs per day
  sta_daily_df<-data.frame(staID=rep(as.character(j),length(sta_data_daily_xts)),
                           date=as.Date(strptime(index(sta_data_daily_xts),format="%Y-%m-%d %H:%M"),format="%Y-%m-%d"),
                           obs_int_mins=rep(obs_int_minutes,length(sta_data_daily_xts)),
                           data_per=sta_per_obs_daily_xts,
                           rf=sta_data_daily_xts) #make df row
  meso_daily_rf<-rbind(meso_daily_rf,sta_daily_df)
}
print("loop complete!")
row.names(meso_daily_rf)<-NULL #rename rows

#head(meso_daily_rf)
#tail(meso_daily_rf)

#subsets: yesterday with 95% data
meso_daily_rf_today<-meso_daily_rf[meso_daily_rf$date==(currentDate),]#subset yesterday
meso_daily_rf_today_final<-meso_daily_rf_today[meso_daily_rf_today$data_per>=0.95,]#subset days with at least 95% data
row.names(meso_daily_rf_today_final)<-NULL #rename rows
#head(meso_daily_rf_today)
#tail(meso_daily_rf_today)

#write or append daily rf data monthly file
setwd(agg_daily_wd)#server path daily agg file
rf_month_filename<-paste0(format((currentDate),"%Y_%m"),"_synoMeso_daily_rf.csv") #dynamic file name that includes month year so when month is done new file is written

if(file.exists(rf_month_filename)){
  write.table(meso_daily_rf_today_final,rf_month_filename, row.names=F,sep = ",", col.names = F, append = T)
  print(paste(rf_month_filename,"appended"))
}else{
  write.csv(meso_daily_rf_today_final,rf_month_filename, row.names=F)
  print(paste(rf_month_filename,"written"))
}

print("PAU!")

