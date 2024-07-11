#this code grabs all raw madis data for HI and calcs daily RF total from precip accumulations

rm(list = ls())#remove all objects in R

#set options
options(warn=-1)#supress warnings for session
Sys.setenv(TZ='Pacific/Honolulu') #set TZ to honolulu
print(paste("madis rf daily run:",Sys.time()))#for cron log

#set MAIN DIR
mainDir <- "/home/hawaii_climate_products_container/preliminary"
codeDir<-paste0(mainDir,"/rainfall/code/source")

#define dates
source(paste0(codeDir,"/dataDateFunc.R"))
dataDate<-dataDateMkr() #function for importing/defining date as input or as yesterday
currentDate<-dataDate #dataDate as currentDate (yesterday)

#load packages
#install.packages("xts","lubridate")
require(xts)
require(lubridate)

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
PC2rf24hr<-function(df,pc_dp,doi){
  require(lubridate)
  dtoi1<-strptime(paste(doi, "00:00:00"),format="%Y-%m-%d %H:%M:%S")
  dtoi2<-strptime(paste(doi+1, "00:00:00"),format="%Y-%m-%d %H:%M:%S")
  attr(dtoi1,"tzone") <- "Pacific/Honolulu" #convert TZ attribute to HST
  attr(dtoi2,"tzone") <- "Pacific/Honolulu" #convert TZ attribute to HST
  
  #data wrangle
  df_rp<-subset(df, varname=="rawPrecip")# subset rawPrecip only
  df_rp<-df_rp[complete.cases(df_rp),] #remove NA rows
  df_rp$time<-strptime(df_rp$time, format="%Y-%m-%d %H:%M", tz="UTC")
  attr(df_rp$time,"tzone") <- "Pacific/Honolulu" #convert TZ attribute to HST
  df_rp$time<-(df_rp$time-36000) #minus 10 hrs for UTC to HST
  df_rp<-df_rp[df_rp$dataProvider %in% pc_dp,] #subset data providers with PC var as rainfall measurement
  df_rp$time<-lubridate::round_date(df_rp$time, "10 minutes")
  max(df_rp$time)
  
  #pc1
  df_pc1<-df_rp[df_rp$time==dtoi1,]  #subset only matching date time pc obs
  df_pc1<-df_pc1[,c("stationId","value")]
  names(df_pc1)[ncol(df_pc1)]<-"PC1"
  df_pc1<-df_pc1[order(df_pc1$PC1),] #reorder PC1 smallest to biggest
  df_pc1<-df_pc1[!duplicated(df_pc1$stationId),] #remove dp station with larger PC1
  row.names(df_pc1)<-NULL
  
  #pc2
  df_pc2<-df_rp[df_rp$time==dtoi2,]  #subset only matching date time pc obs
  df_pc2<-df_pc2[,c("stationId","value")]
  names(df_pc2)[ncol(df_pc2)]<-"PC2"
  df_pc2<-df_pc2[order(df_pc2$PC2,decreasing = T),] #reorder PC1  biggest to smallest
  df_pc2<-df_pc2[!duplicated(df_pc2$stationId),] #remove dp station with larger PC1
  row.names(df_pc2)<-NULL
  
  #merge
  mergedPC<-merge(df_pc1,df_pc2,by="stationId")
  mergedPC$PCdiff<-mergedPC$PC2-mergedPC$PC1
  mergedPC<-mergedPC[mergedPC$PCdiff>=0,] #subset minus values
  
  #make matching DF
  sta_daily_df<-data.frame(staID=mergedPC$stationId,
                           date=rep(doi,nrow(mergedPC)),
                           obs_int_mins=NA,
                           data_per=NA,
                           rf=round(mergedPC$PCdiff,4))#make df row
  
  return(sta_daily_df)
}#end PC2rf24hr

#dirs
parse_wd<-paste0(mainDir,"/data_aqs/data_outputs/madis/parse")
agg_daily_wd<-paste0(mainDir,"/rainfall/working_data/madis")

#read MADIS parsed table from dev server OLD
# setwd(parse_wd)#sever path for parsed madis files
# madis_filename<-paste0(format((currentDate),"%Y%m%d"),"_madis_parsed.csv") #dynamic file name that includes date
# all_madis<-read.csv(madis_filename)

#read MADIS parsed table from ikewai data portal
ikeUrl<-"https://ikeauth.its.hawaii.edu/files/v2/download/public/system/ikewai-annotated-data/HCDP/workflow_data/preliminary" #url
madis_filename<-paste0(ikeUrl,"/data_aqs/data_outputs/madis/parse/",format((currentDate),"%Y%m%d"),"_madis_parsed.csv") #dynamic file name that includes date
madis_filename1<-paste0(ikeUrl,"/data_aqs/data_outputs/madis/parse/",format((currentDate-1),"%Y%m%d"),"_madis_parsed.csv") #dynamic file name that includes date-1

all_madis<-rbind(read.csv(madis_filename),read.csv(madis_filename1))#read data date and day before data date csv
#head(all_madis)

#subset precip var, convert inch to mm and convert UTC to HST
all_madis_rp<-subset(all_madis, varname=="rawPrecip")# subset rawPrecip only
all_madis_rp<-all_madis_rp[complete.cases(all_madis_rp),] #remove NA rows
all_madis_rp$time<-strptime(all_madis_rp$time, format="%Y-%m-%d %H:%M", tz="UTC")
attr(all_madis_rp$time,"tzone") <- "Pacific/Honolulu" #convert TZ attribute to HST
all_madis_rp$time<-(all_madis_rp$time-36000) #minus 10 hrs for UTC to HST
all_madis_rp$time<-(all_madis_rp$time)-1 #minus 1 second to put midnight obs in last day
#tail(all_madis_rp)
#head(all_madis_rp)

##incremental rainfall
#subset by data providers with PP (incremental rainfall) in rawPrecip var
pp_dp<-c("HF-METAR")
all_madis_pp<-all_madis_rp[all_madis_rp$dataProvider %in% pp_dp,] #subset data providers with PP var as rainfall measurement

#unique madis stations
stations<-unique(all_madis_pp$stationId)

#blank DF to store daily data
madis_daily_pp_rf<-data.frame()

#start daily RF loop for PP
print("PP daily rf loop started...")
for(j in stations){
  sta_data<-subset(all_madis_pp,stationId==j)
  sta_data_xts<-xts(sta_data$value,order.by=sta_data$time,unique = TRUE) #make xtended timeseries object
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
    sta_daily_df<-data.frame(staID=rep(as.character(j),length(sta_data_daily_xts)),date=as.Date(strptime(index(sta_data_daily_xts),format="%Y-%m-%d %H:%M"),format="%Y-%m-%d"),obs_int_mins=rep(obs_int_minutes,length(sta_data_daily_xts)),data_per=sta_per_obs_daily_xts,rf=sta_data_daily_xts)#make df row
    sta_daily_df$rf<-round(sta_daily_df$rf,4) #round rainfall 
    madis_daily_pp_rf<-rbind(madis_daily_pp_rf,sta_daily_df)
  }
}
print("PP loop complete!")

#head(madis_daily_pp_rf)
#tail(madis_daily_pp_rf)

#subset yesterday 
madis_daily_pp_rf<-madis_daily_pp_rf[madis_daily_pp_rf$date==(currentDate),]#subset yesterday

#subset 95% data
madis_daily_pp_rf<-madis_daily_pp_rf[madis_daily_pp_rf$data_per>=0.95,]#subset days with at least 95% data
madis_daily_pp_rf<-madis_daily_pp_rf[order(madis_daily_pp_rf$data_per,decreasing = T),] #sort descending by data percent
row.names(madis_daily_pp_rf)<-NULL #rename rows

##accumulated rainfall
#subset by data providers with PC (accumulated rainfall) in rawPrecip var
pc_dp<-c("HADS","RAWS","MesoWest")
all_madis_pc<-all_madis_rp[all_madis_rp$dataProvider %in% pc_dp,] #subset data providers with PC var as rainfall measurement

#unique madis stations
stations<-unique(all_madis_pc$stationId)

#blank DF to store daily data
madis_daily_pc_rf<-data.frame()

#start PC to incremental daily RF loop
print("PC daily rf loop started...")
for(j in stations){
  sta_data<-subset(all_madis_pc,stationId==j)
  sta_data_xts<-xts(sta_data$value,order.by=sta_data$time,unique = TRUE) #make xtended timeseries object
  sta_data_xts_sub<- sta_data_xts[!duplicated(index(sta_data_xts)),] #remove duplicate time obs
  if(nrow(sta_data_xts_sub)>=23){ #only stations with at least hourly intervals
    sta_data_xts_sub_lag<-diff(sta_data_xts_sub,lag=1)
    sta_data_xts_sub_lag[sta_data_xts_sub_lag<0]<-NA #NA to neg values when lag 1 dif
    sta_data_hrly_xts<-apply.hourly(sta_data_xts_sub_lag,FUN=sum,roundtime = "trunc")#agg to hourly and truncate hour
    indexTZ(sta_data_hrly_xts) <- "Pacific/Honolulu"
    sta_data_daily_xts<-apply.daily(sta_data_hrly_xts,FUN=sum,na.rm = F)#daily sum of all hrly observations
    obs_ints<-diff(index(sta_data_xts_sub),lag=1) #calculate vector of obs intervals
    obs_int_hr<-getmode(as.numeric(obs_ints, units="hours"))
    obs_int_minutes<-obs_int_hr*60
    obs_per_day<-((1/obs_int_hr)*24)#calculate numbers of obs per day based on obs interval
    sta_per_obs_daily_xts<-as.numeric(apply.daily(sta_data_xts_sub_lag,FUN=length)/obs_per_day)#vec of % percentage of obs per day
    sta_daily_df<-data.frame(staID=rep(as.character(j),length(sta_data_daily_xts)),date=as.Date(strptime(index(sta_data_daily_xts),format="%Y-%m-%d %H:%M"),format="%Y-%m-%d"),obs_int_mins=rep(obs_int_minutes,length(sta_data_daily_xts)),data_per=sta_per_obs_daily_xts,rf=sta_data_daily_xts)#make df row
    sta_daily_df$rf<-round(sta_daily_df$rf,4) #round rainfall 
    madis_daily_pc_rf<-rbind(madis_daily_pc_rf,sta_daily_df)
  }
}
print("PC loop complete!")
#head(madis_daily_pc_rf)
#tail(madis_daily_pc_rf)

#subset yesterday
madis_daily_pc_rf<-madis_daily_pc_rf[madis_daily_pc_rf$date==(currentDate),]#subset data day
#head(madis_daily_pc_rf)
#tail(madis_daily_pc_rf)

#subset 95% data for day
madis_daily_pc_rf<-madis_daily_pc_rf[madis_daily_pc_rf$data_per>=0.95,]#subset days with at least 95% data
row.names(madis_daily_pc_rf)<-NULL #rename rows

#data check
print("data check head:tail..")
head(madis_daily_pc_rf)
tail(madis_daily_pc_rf)

##alternate 24hr PC processing
madis_daily_pc_rf_alt<-PC2rf24hr(df=all_madis,pc_dp=pc_dp,doi=currentDate)

#combining all data and remove duplicates
all_madis_daily_pc<-rbind(madis_daily_pc_rf,madis_daily_pc_rf_alt)#rbind all processed PC data with with alt 24hr pc
all_madis_daily_rf<-rbind(madis_daily_pp_rf,all_madis_daily_pc) #add pp station data
all_madis_daily_rf_final<-all_madis_daily_rf[!duplicated(all_madis_daily_rf$staID),] #remove dup stations
row.names(all_madis_daily_rf_final)<-NULL #rename rows

#final data check
print("final check of data: head,tail...")
head(all_madis_daily_rf_final)
tail(all_madis_daily_rf_final)

#write or append daily rf data monthly file
setwd(agg_daily_wd)#server path daily agg file
rf_month_filename<-paste0(format((currentDate),"%Y_%m"),"_madis_daily_rf.csv") #dynamic file name that includes month year so when month is done new file is written

#local write/append
if(file.exists(rf_month_filename)){
	 write.table(all_madis_daily_rf_final,rf_month_filename, row.names=F,sep = ",", col.names = F, append = T)
	 print(paste(rf_month_filename,"appended"))
      }else{
	 write.csv(all_madis_daily_rf_final,rf_month_filename, row.names=F)
	 print(paste(rf_month_filename,"written"))
      }


print("PAU!")
