#this code grabs all raw HADS data for HI and calcs daily RF total from precip accumulations

rm(list = ls())#remove all objects in R

#set global options
options(warn=-1) #suppress warnings for session
Sys.setenv(TZ='Pacific/Honolulu') #set TZ to honolulu
print(paste("hads rf daily run:",Sys.time()))#for cron log

#load packages
#install.packages("xts")
require(xts)

#set dirs
mainDir <- "/home/hawaii_climate_products_container/preliminary"
codeDir<-paste0(mainDir,"/rainfall/code/source")
parse_wd<-paste0(mainDir,"/data_aqs/data_outputs/hads/parse")
agg_daily_wd<-paste0(mainDir,"/rainfall/working_data/hads")

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
PC2rf24hr<-function(data,doi){
  dtoi1<-strptime(paste(doi, "00:00:00"),format="%Y-%m-%d %H:%M:%S")
  dtoi2<-strptime(paste(doi+1, "00:00:00"),format="%Y-%m-%d %H:%M:%S")
  attr(dtoi1,"tzone") <- "Pacific/Honolulu" #convert TZ attribute to HST
  attr(dtoi2,"tzone") <- "Pacific/Honolulu" #convert TZ attribute to HST
  
  #hads source
    hads_pc<-subset(data,var=="PC")# subset precip only
    hads_pc$random<-trimws(as.character(hads_pc$random))
    hads_pc<-subset(hads_pc, random == "" | is.na(random))[,-c(6)] #remove random samples and random col
    hads_pc$value<-hads_pc$value*25.4 #convert to MM
    hads_pc<-hads_pc[complete.cases(hads_pc),] #remove NA rows
    hads_pc$obs_time<-strptime(hads_pc$obs_time, format="%Y-%m-%d %H:%M", tz="UTC")
    attr(hads_pc$obs_time,"tzone") <- "Pacific/Honolulu" #convert TZ attribute to HST
    hads_pc$obs_time<-(hads_pc$obs_time-36000) #minus 10 hrs for UTC to HST
    #pc1
    hads_pc1<-hads_pc[hads_pc$obs_time==dtoi1,]  #subset only matching date time pc obs
    hads_pc1<-hads_pc1[,c("staID","NWS_sid","value")]
    names(hads_pc1)[ncol(hads_pc1)]<-"PC1"
    row.names(hads_pc1)<-NULL
    #pc2
    hads_pc2<-hads_pc[hads_pc$obs_time==dtoi2,]  #subset only matching date time pc obs
    hads_pc2<-hads_pc2[,c("staID","value")]
    names(hads_pc2)[ncol(hads_pc2)]<-"PC2"
    row.names(hads_pc2)<-NULL
    
    #merge
    mergedPC<-merge(hads_pc1,hads_pc2,by="staID")
    mergedPC$PCdiff<-mergedPC$PC2-mergedPC$PC1
    mergedPC<-mergedPC[mergedPC$PCdiff>=0,] #subset minus values
    
    #make matching DF
    daily_df<-data.frame(NWS_sid=mergedPC$NWS_sid,
                             staID=mergedPC$staID,
                             date=rep(doi,nrow(mergedPC)),
                             obs_int_mins=rep(NA,nrow(mergedPC)),
                             data_per=rep(NA,nrow(mergedPC)),
                             rf=mergedPC$PCdiff)#make df row
    return(daily_df)
}#end PC2rf24hr

# #read HADS parsed table from dev server OLD
# setwd(parse_wd)#sever path for parsed hads files
# hads_filename<-paste0(format((currentDate),"%Y%m%d"),"_hads_parsed.csv") #dynamic file name that includes date
# all_hads<-read.csv(hads_filename)
# #head(all_hads)

#read HADS parsed table from ikewai data portal
ikeUrl<-"https://ikeauth.its.hawaii.edu/files/v2/download/public/system/ikewai-annotated-data/HCDP/workflow_data/preliminary" #url
hads_filename<-paste0(format((currentDate),"%Y%m%d"),"_hads_parsed.csv") #dynamic file name that includes date
all_hads<-read.csv(paste0(ikeUrl,"/data_aqs/data_outputs/hads/parse/",hads_filename))
#head(all_hads)

#subset precip var, convert inch to mm and convert UTC to HST
all_hads_pc_pp<-subset(all_hads,var=="PC"|var=="PP")# subset precip only
all_hads_pc_pp$random<-trimws(as.character(all_hads_pc_pp$random))
all_hads_pc_pp<-subset(all_hads_pc_pp, random == "" | is.na(random))[,-c(6)] #remove random samples and random col
all_hads_pc_pp$value<-all_hads_pc_pp$value*25.4 #convert to MM
all_hads_pc_pp<-all_hads_pc_pp[complete.cases(all_hads_pc_pp),] #remove NA rows
all_hads_pc_pp$obs_time<-strptime(all_hads_pc_pp$obs_time, format="%Y-%m-%d %H:%M", tz="UTC")
attr(all_hads_pc_pp$obs_time,"tzone") <- "Pacific/Honolulu" #convert TZ attribute to HST
all_hads_pc_pp$obs_time<-(all_hads_pc_pp$obs_time-36000) #minus 10 hrs for UTC to HST
all_hads_pc_pp$obs_time<-(all_hads_pc_pp$obs_time)-1 #minus 1 second to put midnight obs in last day
#tail(all_hads_pc_pp,3)
#head(all_hads_pc_pp,3)

## process PP vars
all_hads_pp<-subset(all_hads_pc_pp,var=="PP")# subset precip only
str(all_hads_pp)

#blank DF to store daily data
hads_daily_pp_rf<-data.frame()

#unique hads stations
stations<-unique(all_hads_pp$staID)

#start daily RF loop for PP
print("daily rf loop started...")
for(j in stations){
  sta_data<-subset(all_hads_pp,staID==j)
  sta_data_xts<-xts(sta_data$value,order.by=sta_data$obs_time,unique = TRUE) #make xtended timeseries object
  sta_data_xts_sub<- sta_data_xts[!duplicated(index(sta_data_xts)),] #remove duplicate time obs
  if(nrow(sta_data_xts_sub)>=23){
    sta_data_hrly_xts<-apply.hourly(sta_data_xts_sub,FUN=sum,roundtime = "trunc")#agg to hourly and truncate hour
    indexTZ(sta_data_hrly_xts) <- "Pacific/Honolulu"
    sta_data_daily_xts<-apply.daily(sta_data_hrly_xts,FUN=sum,na.rm = F)#daily sum of all all lag observations
    obs_ints<-diff(index(sta_data_xts_sub),lag=1) #calculate vector of obs intervals
    obs_int_hr<-getmode(as.numeric(obs_ints, units="hours"))
    obs_int_minutes<-obs_int_hr*60
    obs_per_day<-((1/obs_int_hr)*24)#calculate numbers of obs per day based on obs interval
    sta_per_obs_daily_xts<-as.numeric(apply.daily(sta_data_xts_sub,FUN=length)/obs_per_day)#vec of % percentage of obs per day
    sta_daily_df<-data.frame(NWS_sid=rep(as.character(unique(sta_data$NWS_sid)),length(sta_data_daily_xts)),staID=rep(as.character(j),length(sta_data_daily_xts)),date=as.Date(strptime(index(sta_data_daily_xts),format="%Y-%m-%d %H:%M"),format="%Y-%m-%d"),obs_int_mins=rep(obs_int_minutes,length(sta_data_daily_xts)),data_per=sta_per_obs_daily_xts,rf=sta_data_daily_xts)#make df row
    hads_daily_pp_rf<-rbind(hads_daily_pp_rf,sta_daily_df)
  }
}
print("loop complete!")

#head(hads_daily_pp_rf)
#tail(hads_daily_pp_rf)

#subset yesterday 
hads_daily_pp_rf<-hads_daily_pp_rf[hads_daily_pp_rf$date==(currentDate),]#subset yesterday

#subset 95% data
hads_daily_pp_rf<-hads_daily_pp_rf[hads_daily_pp_rf$data_per>=0.95,]#subset days with at least 95% data
hads_daily_pp_rf<-hads_daily_pp_rf[order(hads_daily_pp_rf$data_per,decreasing = T),] #sort descending by data percent
row.names(hads_daily_pp_rf)<-NULL #rename rows


##process PC vars
all_hads_pc<-subset(all_hads_pc_pp,var=="PC")# subset precip only
str(all_hads_pc)
all_hads_pc[all_hads_pc$staID=="1560740C","obs_time"]

#blank DF to store daily data
hads_daily_pc_rf<-data.frame()

#unique hads stations
stations<-unique(all_hads_pc$staID)

#start daily RF loop for PC
print("daily rf loop started...")
for(j in stations){
  sta_data<-subset(all_hads_pc,staID==j)
  sta_data_xts<-xts(sta_data$value,order.by=sta_data$obs_time,unique = TRUE) #make xtended timeseries object
  sta_data_xts_sub<- sta_data_xts[!duplicated(index(sta_data_xts)),] #remove duplicate time obs
  if(nrow(sta_data_xts_sub)>=23){
    sta_data_xts_sub_lag<-diff(sta_data_xts_sub,lag=1)
    sta_data_xts_sub_lag[sta_data_xts_sub_lag<0]<-NA #NA to neg values when lag 1 dif
    sta_data_hrly_xts<-apply.hourly(sta_data_xts_sub_lag,FUN=sum,roundtime = "trunc")#agg to hourly and truncate hour
    indexTZ(sta_data_hrly_xts) <- "Pacific/Honolulu"
    sta_data_daily_xts<-apply.daily(sta_data_hrly_xts,FUN=sum,na.rm = F)#daily sum of all all lag observations
    obs_ints<-diff(index(sta_data_xts_sub),lag=1) #calculate vector of obs intervals
    obs_int_hr<-getmode(as.numeric(obs_ints, units="hours"))
    obs_int_minutes<-obs_int_hr*60
    obs_per_day<-((1/obs_int_hr)*24)#calculate numbers of obs per day based on obs interval
    sta_per_obs_daily_xts<-as.numeric(apply.daily(sta_data_xts_sub_lag,FUN=length)/obs_per_day)#vec of % percentage of obs per day
    sta_daily_df<-data.frame(NWS_sid=rep(as.character(unique(sta_data$NWS_sid)),length(sta_data_daily_xts)),staID=rep(as.character(j),length(sta_data_daily_xts)),date=as.Date(strptime(index(sta_data_daily_xts),format="%Y-%m-%d %H:%M"),format="%Y-%m-%d"),obs_int_mins=rep(obs_int_minutes,length(sta_data_daily_xts)),data_per=sta_per_obs_daily_xts,rf=sta_data_daily_xts)#make df row
    hads_daily_pc_rf<-rbind(hads_daily_pc_rf,sta_daily_df)
  }
}
print("loop complete!")

#head(hads_daily_pc_rf)
#tail(hads_daily_pc_rf)

#subset yesterday 
hads_daily_pc_rf<-hads_daily_pc_rf[hads_daily_pc_rf$date==(currentDate),]#subset yesterday

#subset 95% data
hads_daily_pc_rf<-hads_daily_pc_rf[hads_daily_pc_rf$data_per>=0.95,]#subset days with at least 95% data
hads_daily_pc_rf<-hads_daily_pc_rf[order(hads_daily_pc_rf$data_per,decreasing = T),] #sort descending by data percent
row.names(hads_daily_pc_rf)<-NULL #rename rows

##alternate 24hr PC processing
hads_daily_pc_rf_alt<-PC2rf24hr(data=all_hads,doi=currentDate)

all_hads_daily_pc<-rbind(rbind(hads_daily_pc_rf[hads_daily_pc_rf$data_per==1,],hads_daily_pc_rf_alt),hads_daily_pc_rf[hads_daily_pc_rf$data_per<1,]) #rbind all processed PC data with 100% data then, alt 24hr pc, then pc with 95%> & <100%
all_hads_daily_rf<-rbind(all_hads_daily_pc,hads_daily_pp_rf) #add pp station data
all_hads_daily_rf_final<-all_hads_daily_rf[!duplicated(all_hads_daily_rf$staID),] #remove dup stations
row.names(all_hads_daily_rf_final)<-NULL #rename rows

#final data check
#print(all_hads_daily_rf_final)
head(all_hads_daily_rf_final)
tail(all_hads_daily_rf_final)

#write or append daily rf data monthly file
setwd(agg_daily_wd)#server path daily agg file
rf_month_filename<-paste0(format((currentDate),"%Y_%m"),"_hads_daily_rf.csv") #dynamic file name that includes month year so when month is done new file is written

if(file.exists(rf_month_filename)){
  write.table(all_hads_daily_rf_final,rf_month_filename, row.names=F,sep = ",", col.names = F, append = T)
  print(paste(rf_month_filename,"appended"))
}else{
  write.csv(all_hads_daily_rf_final,rf_month_filename, row.names=F)
  print(paste(rf_month_filename,"written"))
}

print("PAU!")