##check daily rf station qa/qc and give prob/flag "bad" values
#install.packages(c("reshape2","raster","rgdal","randomForest","geosphere","data.table"))
rm(list=ls())

library(reshape2)
library(raster)
library(rgdal)
library(geosphere)
library(data.table)
library(randomForest)

#set MAIN DIR
mainDir <- "/home/hawaii_climate_products_container/preliminary"
codeDir<-paste0(mainDir,"/rainfall/code/source")

#input dirs
rf_day_data_wd <- paste0(mainDir,"/rainfall/data_outputs/tables/station_data/daily/raw/statewide") #input final combine daily rainfall data
lognorm_ras_wd<- paste0(mainDir,"/rainfall/dependencies/daily/probRasters") #input lognorm ras depend
qaqc_models_wd<- paste0(mainDir,"/rainfall/dependencies/daily/models") #input mod depend

#output dirs
qaqc_fail_stations_wd<- paste0(mainDir,"/rainfall/data_outputs/tables/rf_station_tracking/qaqc_fail")
rf_qaqc_flag_data_wd<- paste0(mainDir,"/rainfall/data_outputs/tables/station_data/daily/qc_flag/statewide")
rf_qaqc_prob_data_wd<- paste0(mainDir,"/rainfall/data_outputs/tables/station_data/daily/qc_prob/statewide")
rf_qaqc_screen_data_wd<- paste0(mainDir,"/rainfall/data_outputs/tables/station_data/daily/raw_qc/statewide")
rf_qaqc_count_data_wd <- paste0(mainDir,"/rainfall/data_outputs/tables/rf_station_tracking/qaqc_count")

#rf to prob function
rf_Prob<-function(rf,meanlog,sdlog,pop){
  if(rf>0){
    rf_per<-plnorm(rf,meanlog,sdlog, lower.tail = T)#get lower tail of daily rf val (ie:cumulative prob of rf being this value or lower)
    return((1-pop)+(pop*rf_per))#include 1-pop into prob
  }else{ #if rf is zero
    return(1-pop) #assign 1-pop ie: prob of it not raining 
  }#end if/else
}#end function

#calc log prob
get_LogProb<-function(sta_df,lnAnnRas,ln6moRas){
  coordinates(sta_df) = ~LON + LAT #make sta df spatial
  log_norm_sta<-as.data.frame(sta_df,extract(lnAnnRas,sta_df),extract(ln6moRas,sta_df))
  names(log_norm_sta)[4]<-"rf_mm"
  #calc ann & 6 mon rf prob
  sta_df$rf_6mo_per<-mapply(rf_Prob,rf=log_norm_sta$rf_mm,meanlog=log_norm_sta$meanlog_6mo,sdlog=log_norm_sta$sdlog_6mo, pop=log_norm_sta$pop_6mo)# convert random rf values	to 6mo per	
  sta_df$rf_ann_per<-mapply(rf_Prob,rf=log_norm_sta$rf_mm,meanlog=log_norm_sta$meanlog_ann,sdlog=log_norm_sta$sdlog_ann, pop=log_norm_sta$pop_ann)# convert random rf values	to ann per
  final_df<-as.data.frame(sta_df)[,c(1,4:6)]
  names(final_df)<-c("SKN_target","rf_target","rf_6mo_per_target","rf_ann_per_target")
  return(final_df)
}

#calc close 10 with probs
close10RF<-function(day_df,ln6moRas,lnAnnRas){
  dist <- as.data.frame(distm(day_df[,c("LON","LAT")]))#dist of all stations to all other stations for day col 'i'
  names(dist)<-as.character(day_df$SKN)
  #row.names(dist)<-as.character(day_df$SKN)
  dist$SKN<-as.character(day_df$SKN)
  dist_long<-reshape2::melt(dist, id.vars=c("SKN"))
  names(dist_long)<-c("SKN_obs","SKN_target","dist")
  dist_long$SKN_obs<-as.character(dist_long$SKN_obs)
  dist_long$SKN_target<-as.character(dist_long$SKN_target)
  targets<-unique(dist_long$SKN_target)
  dist10_day<-data.frame()
  for(p in targets){
    dist_SKN<-dist_long[dist_long$SKN_target==p,]
    dist_SKN<-dist_SKN[dist_SKN$SKN_target!=dist_SKN$SKN_obs,]
    dist10_SKN<-dist_SKN[(order(dist_SKN[,"dist"])[1:10]),]#10 closest stations
    row.names(dist10_SKN)<-NULL
    dist10_SKN$dist_rank<-row.names(dist10_SKN)
    dist10_SKN<-dist10_SKN[,c(2,1,3,4)]
    dist10_SKN$date<-rep(as.character(names(day_df)[4]))
    dist10_SKN$date_skn_target<-paste(dist10_SKN$date,dist10_SKN$SKN_target,sep="_")
    dist10_SKN$date_skn_obs<-paste(dist10_SKN$date,dist10_SKN$SKN_obs,sep="_")
    dist10_SKN$date_skn_tar_obs<-paste(dist10_SKN$date,dist10_SKN$SKN_obs,dist10_SKN$SKN_target,sep="_")
    dist10_day<-as.data.frame(rbind(dist10_day,dist10_SKN))
  }
  
  day_df_prob<-get_LogProb(day_df,lnAnnRas,ln6moRas)
  #str(day_df_prob)
  dist10_day_merge<-merge(dist10_day,day_df_prob,by="SKN_target")
  names(day_df_prob)<-gsub("target","obs",names(day_df_prob))
  dist10_day_long<-merge(dist10_day_merge,day_df_prob,by="SKN_obs")
  dist10_day_wide <- data.table::dcast(as.data.table(dist10_day_long), SKN_target + rf_target + rf_6mo_per_target + rf_ann_per_target ~ dist_rank, value.var = c("dist","rf_obs","rf_6mo_per_obs","rf_ann_per_obs"))
  names(day_df)[1]<-"SKN_target"
  final10dist<-merge(day_df[,c("SKN_target","LON","LAT")],as.data.frame(dist10_day_wide),by="SKN_target")
  final10dist<-final10dist[,c("SKN_target","LON","LAT","rf_target","rf_6mo_per_target","rf_ann_per_target","dist_1","dist_2","dist_3","dist_4","dist_5","dist_6","dist_7","dist_8","dist_9","rf_obs_1","rf_obs_10","rf_obs_2","rf_obs_3","rf_obs_4","rf_obs_5","rf_obs_6","rf_obs_7","rf_obs_8","rf_obs_9","dist_10",
                              "rf_6mo_per_obs_1","rf_6mo_per_obs_2","rf_6mo_per_obs_3","rf_6mo_per_obs_4","rf_6mo_per_obs_5","rf_6mo_per_obs_6","rf_6mo_per_obs_7","rf_6mo_per_obs_8","rf_6mo_per_obs_9","rf_6mo_per_obs_10",
                              "rf_ann_per_obs_1","rf_ann_per_obs_2","rf_ann_per_obs_3","rf_ann_per_obs_4","rf_ann_per_obs_5","rf_ann_per_obs_6","rf_ann_per_obs_7","rf_ann_per_obs_8","rf_ann_per_obs_9","rf_ann_per_obs_10")]
  return(final10dist)
}

#bad pred function
pred_bad<-function(close10DF,rfDayCol,rfZ,rfNZ){
  close10DFNZ<-close10DF[close10DF$rf_target>0,]
  close10DFZERO<-close10DF[close10DF$rf_target==0,]
  close10DFNZ$statusPred<-as.numeric(as.character(predict(rfNZ,as.data.frame(close10DFNZ))))
  close10DFNZ$probBad<-as.data.frame(predict(rfNZ,as.data.frame(close10DFNZ),type="prob"))[,"1"]
  close10DFZERO$statusPred<-as.numeric(as.character(predict(rfZ,as.data.frame(close10DFZERO))))
  close10DFZERO$probBad<-as.data.frame(predict(rfZ,as.data.frame(close10DFZERO),type="prob"))[,"1"]
  predDF<-rbind(close10DFNZ[,c("SKN_target","statusPred","probBad")],close10DFZERO[,c("SKN_target","statusPred","probBad")])
  finalPredDF<-merge(close10DF[,c("SKN_target","LON","LAT","rf_target")],predDF,by="SKN_target")
  names(finalPredDF)[c(1,4)]<-c("SKN",rfDayCol)
  return(finalPredDF)
}

#get metadata
meta_url <- "https://raw.githubusercontent.com/ikewai/hawaii_wx_station_mgmt_container/main/Hawaii_Master_Station_Meta.csv"
geog_meta<-read.csv(meta_url, colClasses=c("NESDIS.id"="character"))

#define date
source(paste0(codeDir,"/dataDateFunc.R"))
dataDate<-dataDateMkr() #function for importing/defining date as input or as yesterday
map_date<-dataDate #dataDate as map_date
file_date<-format(map_date,"%Y_%m")

#get month data
setwd(rf_day_data_wd)
rf_month_filename<-paste0("Statewide_Raw_Daily_RF_mm_",file_date,".csv") #dynamic file name that includes month year so when month is done new file is written
monthly_rf<-read.csv(rf_month_filename)
monthly_rf$co<-as.factor(monthly_rf$Island)
levels(monthly_rf$co)[levels(monthly_rf$co)=="MA"|levels(monthly_rf$co)=="KO"|levels(monthly_rf$co)=="LA"|levels(monthly_rf$co)=="MO"]<-"MN"
monthly_rf[,grep("X",names(monthly_rf))][monthly_rf[,grep("X",names(monthly_rf))]>500]<-NA #remove daily vals >700 mm
str(monthly_rf)
head(monthly_rf)

#subset day
rfcol<-format(as.Date(map_date),"X%Y.%m.%d")
vars<-c("SKN","co","LAT","LON",rfcol)
if(!as.logical(max(names(monthly_rf)==rfcol))){print(paste("day:",rfcol,"not present in monthly file:",rf_month_filename))}
daily_rf<-monthly_rf[!is.na(monthly_rf[,rfcol]),vars]
row.names(daily_rf)<-NULL
str(daily_rf)

#get mean log sd log & pop rasters
setwd(lognorm_ras_wd)
state_ann<-stack(list.files("ann",full.names = T))
names(state_ann)<-c("meanlog_ann","pop_ann","sdlog_ann")
is_mjjaso<-(as.numeric(format(as.Date(map_date),"%m"))>=5 & as.numeric(format(as.Date(map_date),"%m"))<=10)
state_6mo<-if(is_mjjaso){stack(Sys.glob("6mo/hi_statewide_*mjjaso.tif"))}else{stack(Sys.glob("6mo/hi_statewide_*ndjfma.tif"))}
names(state_6mo)<-c("meanlog_6mo","pop_6mo","sdlog_6mo")

#make blank qaqc df
daily_rf_qaqc_pass<-data.frame()#blank df to store final data checked
daily_rf_qaqc_fail<-data.frame()#blank df to store removed data

s<-Sys.time()
#flag remove repeat bad data
for(i in unique(daily_rf$co)){
  daily_rf_co<-daily_rf[daily_rf[!is.na(daily_rf[,rfcol]),"co"]==i,-2]
  setwd(qaqc_models_wd)
  modNZ<-readRDS(paste0(i,"_no_human_rf_nonzero_randfor.rds"))
  modZ<-readRDS(paste0(i,"_no_human_rf_zero_randfor.rds"))
  close10DF<-close10RF(day_df=daily_rf_co,ln6moRas=state_6mo,lnAnnRas=state_ann)
  daily_rf_co_pred<-pred_bad(close10DF,rfcol,modZ,modNZ)
  if(daily_rf_co_pred[which.max(daily_rf_co_pred$probBad),"statusPred"]==1){
    daily_rf_qaqc_fail<-rbind(daily_rf_qaqc_fail,daily_rf_co_pred[which.max(daily_rf_co_pred$probBad),])
    daily_rf_co_pred<-daily_rf_co_pred[-which.max(daily_rf_co_pred$probBad),]
  }
  nbad<-sum(daily_rf_co_pred$statusPred)
  while(nbad>0){
    daily_rf_co<-daily_rf_co_pred[,1:4]
    close10DF<-close10RF(day_df=daily_rf_co,ln6moRas=state_6mo,lnAnnRas=state_ann)
    daily_rf_co_pred<-pred_bad(close10DF,rfcol,modZ,modNZ)
    if(daily_rf_co_pred[which.max(daily_rf_co_pred$probBad),"statusPred"]==1){
      daily_rf_qaqc_fail<-rbind(daily_rf_qaqc_fail,daily_rf_co_pred[which.max(daily_rf_co_pred$probBad),])
      daily_rf_co_pred<-daily_rf_co_pred[-which.max(daily_rf_co_pred$probBad),]
    }
    nbad<-sum(daily_rf_co_pred$statusPred)
    if(nrow(daily_rf_co_pred)<=11){nbad<-0}
  }#remove bad while loop
  daily_rf_qaqc_pass<-rbind(daily_rf_qaqc_pass,daily_rf_co_pred)
}#co loop end
e<-Sys.time()
print("rf daily value screening...")
print(e-s)

#make qaqc removed station long df
if(ncol(daily_rf_qaqc_fail)==0){
  #blank dataframe to save
  rf_qaqc_fail_stations<-data.frame(SKN=as.numeric(),
                                    date=as.character(),
                                    rfmm=as.numeric(),
                                    statusPred=as.logical(),
                                    probBad=as.numeric())
  rf_qaqc_fail_stations<-merge(geog_meta,rf_qaqc_fail_stations,by="SKN")
  
  }else{
  #make dataframe of fail/flag stations
  rf_qaqc_fail_stations<-data.frame(SKN=daily_rf_qaqc_fail$SKN,
                              date=as.character(rep(as.Date(gsub("X","",rfcol),format="%Y.%m.%d"),ncol=daily_rf_qaqc_fail)),
                              rfmm=daily_rf_qaqc_fail[,rfcol],
                              statusPred=daily_rf_qaqc_fail$statusPred,
                              probBad=daily_rf_qaqc_fail$probBad)

  rf_qaqc_fail_stations<-merge(geog_meta,rf_qaqc_fail_stations,by="SKN")
}
#save or append to monthly qaqc station removes stations file
setwd(qaqc_fail_stations_wd)
rf_fail_stations_month_filename<-paste0(file_date,"_qaqc_fail_daily_rf.csv") #dynamic file name that includes month year so when month is done new file is written

#conditional statement that adds new qaqc daily data
if(file.exists(rf_fail_stations_month_filename)){
  rf_month_fail_df<-read.csv(rf_fail_stations_month_filename,colClasses=c("NESDIS.id"="character"))
  rf_month_fail_df<-rbind(rf_month_fail_df,rf_qaqc_fail_stations)
  write.csv(rf_month_fail_df,rf_fail_stations_month_filename,row.names = F)
  print(paste(rf_fail_stations_month_filename,"daily qaqc fail table appended!"))
}else{
  write.csv(rf_qaqc_fail_stations,rf_fail_stations_month_filename,row.names = F)
  print(paste(rf_fail_stations_month_filename,"daily qaqc fail table written!"))
}

#combine all qaqc results together
daily_rf_checked<-rbind(daily_rf_qaqc_pass,daily_rf_qaqc_fail) #passed and failed stations
daily_rf_checked<-daily_rf_checked[order(as.numeric(daily_rf_checked$SKN)),]
row.names(daily_rf_checked)<-NULL
# str(daily_rf_checked)
# head(daily_rf_checked)

##flag file write or append daily rf QAQC
setwd(rf_qaqc_flag_data_wd)
rf_qaqc_flag_month_filename<-paste0("Statewide_RF_QAQC_Daily_Bad_Flag_",file_date,".csv") #dynamic file name that includes month year so when month is done new file is written

#conditional statement that adds new qaqc daily data
if(file.exists(rf_qaqc_flag_month_filename)){
  #load flag day file add all metadata
  rf_month_flag_df<-read.csv(rf_qaqc_flag_month_filename)
  sub_cols<-c("SKN",names(rf_month_flag_df)[grep("X",names(rf_month_flag_df))])
  rf_month_flag_df<-rf_month_flag_df[,sub_cols]
  rf_month_flag_df_merged<-merge(geog_meta,rf_month_flag_df,by="SKN",all=T)
  #add day flag col to file
  daily_rf_flag<-daily_rf_checked[,c("SKN","statusPred")] #get skn and flag cols
  names(daily_rf_flag)[2]<-rfcol #rename flag col as date
  final_rf_daily_month_qaqc_flag<-merge(rf_month_flag_df_merged,daily_rf_flag,by="SKN",all=T)
  dateCols<-grep("X",names(final_rf_daily_month_qaqc_flag))
  allNA <- apply(final_rf_daily_month_qaqc_flag[,dateCols], 1, function(x) all(is.na(x)))
  final_rf_daily_month_qaqc_flag<-final_rf_daily_month_qaqc_flag[!allNA,]
  write.csv(final_rf_daily_month_qaqc_flag,rf_qaqc_flag_month_filename, row.names=F)
  print(paste(rf_qaqc_flag_month_filename,"daily rf flag table appended!"))
  }else{ #if month year file does not exist make a new month year file
  #new flag file
  daily_rf_flag<-daily_rf_checked[,c("SKN","statusPred")]
  names(daily_rf_flag)[2]<-rfcol
  final_rf_daily_month_qaqc_flag<-merge(geog_meta,daily_rf_flag,by="SKN")
  write.csv(final_rf_daily_month_qaqc_flag,rf_qaqc_flag_month_filename, row.names=F)
  print(paste(rf_qaqc_flag_month_filename,"daily rf flag table written!"))
}

##prob file write or append daily rf QAQC
setwd(rf_qaqc_prob_data_wd) 
rf_qaqc_prob_month_filename<-paste0("Statewide_RF_QAQC_Daily_Bad_Prob_",file_date,".csv") #dynamic file name that includes month year so when month is done new file is written

#conditional statement that adds new qaqc daily data
if(file.exists(rf_qaqc_prob_month_filename)){
  #add prob bad day col to file
  rf_month_prob_df<-read.csv(rf_qaqc_prob_month_filename)
  sub_cols<-c("SKN",names(rf_month_prob_df)[grep("X",names(rf_month_prob_df))])
  rf_month_prob_df<-rf_month_prob_df[,sub_cols]
  rf_month_prob_df_merged<-merge(geog_meta,rf_month_prob_df,by="SKN",all=T)
  #add day prob
  daily_rf_prob<-daily_rf_checked[,c("SKN","probBad")] #get skn and prob cols
  names(daily_rf_prob)[2]<-rfcol #rename prob col as date
  final_rf_daily_month_qaqc_prob<-merge(rf_month_prob_df_merged,daily_rf_prob,by="SKN",all=T)
  dateCols<-grep("X",names(final_rf_daily_month_qaqc_prob))
  allNA <- apply(final_rf_daily_month_qaqc_prob[,dateCols], 1, function(x) all(is.na(x)))
  final_rf_daily_month_qaqc_prob<-final_rf_daily_month_qaqc_prob[!allNA,]
  write.csv(final_rf_daily_month_qaqc_prob,rf_qaqc_prob_month_filename, row.names=F)
  print(paste(rf_qaqc_prob_month_filename,"daily rf prob table appended!"))
  }else{ #if month year file does not exist make a new month year file
  #new prob file
  daily_rf_prob<-daily_rf_checked[,c("SKN","probBad")]
  names(daily_rf_prob)[2]<-rfcol
  final_rf_daily_month_qaqc_prob<-merge(geog_meta,daily_rf_prob,by="SKN")
  write.csv(final_rf_daily_month_qaqc_prob,rf_qaqc_prob_month_filename, row.names=F)
  print(paste(rf_qaqc_prob_month_filename,"daily rf prob table written!"))
}

##save screened station values write/append
setwd(rf_qaqc_screen_data_wd) #set qaqc rainfall output wd
rf_qaqc_screen_month_filename<-paste0("Statewide_QAQC_Raw_Daily_RF_mm_",file_date,".csv") #dynamic file name that includes month year so when month is done new file is written

if(file.exists(rf_qaqc_screen_month_filename)){
  #append qaqc rf
  rf_month_screened<-read.csv(rf_qaqc_screen_month_filename)
  sub_cols<-c("SKN",names(rf_month_screened)[grep("X",names(rf_month_screened))])
  rf_month_screened_merged<-merge(geog_meta,rf_month_screened[,sub_cols],by="SKN",all=T)
  #append screened rf col
  daily_rf_screened<-daily_rf_qaqc_pass[,c("SKN",rfcol)]
  final_rf_daily_month_qaqc_flag<-merge(rf_month_screened_merged,daily_rf_screened,by="SKN",all=T)
  dateCols<-grep("X",names(final_rf_daily_month_qaqc_flag))
  allNA <- apply(final_rf_daily_month_qaqc_flag[,dateCols], 1, function(x) all(is.na(x)))
  final_rf_daily_month_qaqc_flag<-final_rf_daily_month_qaqc_flag[!allNA,]
  write.csv(final_rf_daily_month_qaqc_flag,rf_qaqc_screen_month_filename, row.names=F)
  print(paste(rf_qaqc_screen_month_filename,"daily screened rainfall table appended!"))
}else{ #if month year file does not exist make a new month year file
  #write new qaqc
  daily_rf_screened<-daily_rf_qaqc_pass[,c("SKN",rfcol)]
  final_rf_daily_QAQC<-merge(geog_meta,daily_rf_screened,by="SKN")
  write.csv(final_rf_daily_QAQC,rf_qaqc_screen_month_filename, row.names=F)
  print(paste(rf_qaqc_screen_month_filename,"daily screened rainfall table written!"))
}

#make and save pre/post qaqc station data
rawStaCount<-sum(!is.na(daily_rf[,rfcol]))#save pre screen numbers
qaqcStaCount<-sum(!is.na(daily_rf_screened[,rfcol])) #post qaqc station count
qaStaLogRow<-data.frame(date=map_date,rawStations=rawStaCount,qaqcStations=qaqcStaCount,removeSta=(rawStaCount-qaqcStaCount),perRemove=round((rawStaCount-qaqcStaCount)/rawStaCount,3))

##save station count log write/append
setwd(rf_qaqc_count_data_wd) #set qaqc rainfall output wd
rf_qaqc_log_month_filename<-paste0("Statewide_QAQC_Daily_RF_mm_count_log_",file_date,".csv") #dynamic file name that includes month year so when month is done new file is written
if(file.exists(rf_qaqc_log_month_filename)){
  write.table(qaStaLogRow, file = rf_qaqc_log_month_filename, sep = ",",append=T,quote=F,col.names=F,row.names=F)
  print(paste(rf_qaqc_log_month_filename,"table appended!"))
  }else{
  write.csv(qaStaLogRow,rf_qaqc_log_month_filename,row.names=F)
  print(paste(rf_qaqc_log_month_filename,"table written!"))
  }

#code pau
