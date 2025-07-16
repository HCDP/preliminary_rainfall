#this code combines all daily agg rf together and removes duplicates by priority ranking
rm(list = ls())#remove all objects in R

options(warn=-1)#suppress warnings for session

print(paste("all data daily merge run:",Sys.time()))#for cron log

#set dirs
mainDir <- "/home/hawaii_climate_products_container/preliminary"
codeDir<-paste0(mainDir,"/rainfall/code/source")

#define date
source(paste0(codeDir,"/dataDateFunc.R"))
dataDate<-dataDateMkr() #function for importing/defining date as input or as yesterday
map_date<-dataDate #dataDate as map_date
file_date<-format(map_date,"%Y_%m")

#input dirs
hads_daily_wd <- paste0(mainDir,"/rainfall/working_data/hads") #hads daily agg data wd
nws_daily_wd <- paste0(mainDir,"/rainfall/working_data/nws_rr5") #nws daily agg wd
scan_daily_wd <- paste0(mainDir,"/rainfall/working_data/scan") #scan daily agg wd
madis_daily_wd <- paste0(mainDir,"/rainfall/working_data/madis") #madis daily agg wd
synoptic_daily_wd <- paste0(mainDir,"/rainfall/working_data/hi_mesonet/synoptic") #hi mesonet synoptic daily agg wd
#WORKING... hi mesonet daily agg wd
#WORKING... smart ala wai daily agg wd

#output dirs
missing_sta_wd <- paste0(mainDir,"/rainfall/data_outputs/tables/rf_station_tracking/missing") #unknown stations
count_log_wd <- paste0(mainDir,"/rainfall/data_outputs/tables/rf_station_tracking/count/daily") #station counts per data stream per day
rf_day_source_wd <- paste0(mainDir,"/rainfall/data_outputs/tables/station_data/daily/source/statewide") #datastream source of data
rf_day_data_wd <- paste0(mainDir,"/rainfall/data_outputs/tables/station_data/daily/raw/statewide") #final combine daily rainfall data

#functions
read.csv.TC<-function(file,HADS=FALSE){
  tryCatch({
    if(HADS) {
	    out <- read.csv(file, header=TRUE, colClasses=c("staID"="character"))
    } else {
	    out <- read.csv(file, header=TRUE)
    }
  }, error = function(e) NULL)
}

combine_data <- function(data_filename, new_data, new_date_col, geog_meta) {
  # Read data from file
  file_data <- read.csv(data_filename)
  # Subset data by SKN and date cols
  sub_cols <- c("SKN", names(file_data)[grep("X", names(file_data))])
  # Remove old data for day being processed from source table cols
  sub_cols <- sub_cols[sub_cols != new_date_col]
  file_data <- file_data[,sub_cols]

  # Merge data from file with new data
  merged_data <- merge(file_data, new_data, by="SKN", all=T)

  # Sort columns to ensure dates are properly ordered
  merged_data <- merged_data[,order(colnames(merged_data))]

  # Add geographical metadata
  merged_data <- merge(geog_meta, merged_data, by="SKN")
  return(merged_data)
}

rbind.all.columns <- function(x, y) {     #function to smart rbind
  x.diff <- setdiff(colnames(x), colnames(y))
  y.diff <- setdiff(colnames(y), colnames(x))
  x[, c(as.character(y.diff))] <- NA 
  y[, c(as.character(x.diff))] <- NA 
  return(rbind(x, y))}

#add master metadata with SKN and lat long
meta_url <- "https://raw.githubusercontent.com/ikewai/hawaii_wx_station_mgmt_container/main/Hawaii_Master_Station_Meta.csv"
geog_meta<-read.csv(meta_url, colClasses=c("NESDIS.id"="character"))

#names(geog_meta)
#str(geog_meta)
geog_meta_sub<-geog_meta[,c("SKN","NESDIS.id","SCAN.id","NWS.id","SMART_NODE_RF.id","LAT","LON")]
geog_meta_sub$no_sourceID<-geog_meta_sub$SKN
#head(geog_meta_sub)
print("master meta added!")

#add HADS data
setwd(hads_daily_wd)#set data source wd
hads_month_filename<-paste0(file_date,"_hads_daily_rf.csv")#dynamic file name that includes month year so when month is done new file is written
hads<-read.csv.TC(hads_month_filename,HADS=TRUE)
if(!is.null(hads) && nrow(hads)>0){  #did HADS month file exist? 
 hads<-hads[,c("staID","date","rf")]
 names(hads)<-c("sourceID","date","x") #note 'source_id" IS "NESDIS.id" for hads
 hads$date<-as.Date(hads$date)
 hads<-hads[hads$data_per>=0.95,]#subset days with at least 95% data
 hads_day<-hads[hads$date==map_date,] #date sub
 if(nrow(hads_day)>0){ #if hads_day has rows/data
  hads_day$date<-format(hads_day$date,"%Y.%m.%d")
  hads_day<-hads_day[,c("sourceID","date","x")]
  hads_day$x<-as.numeric(hads_day$x)
  #head(hads_day)
  hads_wide<- reshape(hads_day, idvar = "sourceID", timevar = "date", direction = "wide")
  hads_wide$datastream<-"hads"
  #tail(hads_wide)
  hads_wide_merged_all<-merge(hads_wide,geog_meta_sub[,c("NESDIS.id","SKN")],by.x="sourceID",by.y="NESDIS.id",all.x=T)
  missing_hads<-hads_wide_merged_all[is.na(hads_wide_merged_all$SKN),c("sourceID","datastream")] #missing stations
  hads_wide_merged<-hads_wide_merged_all[!is.na(hads_wide_merged_all$SKN),] #remove missing stations with no SKN
  names(hads_wide_merged)[2]<-gsub("x.","X",names(hads_wide_merged)[2])#make lower case x to uppercase X for continuity 
  hads_wide_merged<-hads_wide_merged[!is.na(hads_wide_merged[,2]),] #remove NA rf day vals
  #tail(hads_wide_merged)
  count_log_hads<-data.frame(datastream=as.character("hads"),station_count=as.numeric(nrow(hads_wide_merged)),unique=as.logical(0)) #log of stations acquired
  print(paste(hads_month_filename,"found!",nrow(hads_wide_merged),"stations added!"))
  }else{ #else if nrow(hads_day) = 0 IE:no day data
   hads_wide_merged<-data.frame(sourceID=as.character(NA),date_day=as.numeric(NA),datastream=as.character("hads"),SKN=as.numeric(NA))
   names(hads_wide_merged)[2]<-format(dataDate,"X%Y.%m.%d")
   missing_hads<-data.frame(sourceID=as.character(NA),datastream=as.character(NA)) #missing stations blank df
   count_log_hads<-data.frame(datastream=as.character("hads"),station_count=as.numeric(0),unique=as.logical(0))#log of stations acquired
   print(paste("NO HADS DATA:",map_date,"!"))
   }}else{ #else if hads month df is missing make a blank df
    hads_wide_merged<-data.frame(sourceID=as.character(NA),date_day=as.numeric(NA),datastream=as.character("hads"),SKN=as.numeric(NA))
    names(hads_wide_merged)[2]<-format(dataDate,"X%Y.%m.%d")
	  missing_hads<-data.frame(sourceID=as.character(NA),datastream=as.character(NA)) #missing stations blank df
	  count_log_hads<-data.frame(datastream=as.character("hads"),station_count=as.numeric(0),unique=as.logical(0))#log of stations acquired
    print(paste(hads_month_filename,"MISSING empty DF made!"))
    }
print("hads pau!")

#add NWS data 
setwd(nws_daily_wd)#set data source wd
nws_month_filename<-paste0(file_date,"_nws_daily_rf.csv")#dynamic file name that includes month year so when month is done new file is written
nws<-read.csv.TC(nws_month_filename)
if(!is.null(nws) && nrow(nws)>0){  #did nws month file exist? TRUE
 nws<-nws[,c("nwsli","date","prec_mm_24hr")]
 names(nws)<-c("sourceID","date","x")
 nws$date<-as.Date(nws$date)#format as date
 nws<-nws[nws$hour_count >= 23,] #subset by stations with at least 23 hourly obs (95% ie: only 1 missing hour)
 nws_day<-nws[nws$date==map_date,]
 if(nrow(nws_day)>0){ #if nws_day has rows/data
  nws_day$date<-format(nws_day$date,"%Y.%m.%d")
  nws_day$x<-as.numeric(nws_day$x)
  nws_wide<- reshape(nws_day, idvar = "sourceID", timevar = "date", direction = "wide")
  nws_wide$datastream<-"nws"
  #head(nws_wide)
  #tail(nws_wide)
  nws_wide_merged_all<-merge(nws_wide,geog_meta_sub[,c("NWS.id","SKN")],by.x="sourceID",by.y="NWS.id",all.x=T)
  missing_nws<-nws_wide_merged_all[is.na(nws_wide_merged_all$SKN),c("sourceID","datastream")] #missing stations
  nws_wide_merged<-nws_wide_merged_all[!is.na(nws_wide_merged_all$SKN),] #remove missing stations with no SKN
  names(nws_wide_merged)[2]<-gsub("x.","X",names(nws_wide_merged)[2])#make lower case x to uppercase X for continuity
  nws_wide_merged<-nws_wide_merged[!is.na(nws_wide_merged[,2]),] #remove NA rf day vals
  #tail(nws_wide_merged)
  count_log_nws<-data.frame(datastream=as.character("nws"),station_count=as.numeric(nrow(nws_wide_merged)),unique=as.logical(0))
  print(paste(nws_month_filename,"found!",nrow(nws_wide_merged),"stations added!"))
  }else{ #else if nrow(nws_day) = 0
   nws_wide_merged<-data.frame(sourceID=as.character(NA),date_day=as.numeric(NA),datastream=as.character("nws"),SKN=as.numeric(NA))
   names(nws_wide_merged)[2]<-format(dataDate,"X%Y.%m.%d")
   missing_nws<-data.frame(sourceID=as.character(NA),datastream=as.character(NA)) #missing stations blank df
   count_log_nws<-data.frame(datastream=as.character("nws"),station_count=as.numeric(0),unique=as.logical(0))
   print(paste("NO NWS DATA:",map_date,"!"))
   }}else{ #else if nws month df is missing make a blank df
    nws_wide_merged<-data.frame(sourceID=as.character(NA),date_day=as.numeric(NA),datastream=as.character("nws"),SKN=as.numeric(NA))
    names(nws_wide_merged)[2]<-format(dataDate,"X%Y.%m.%d")
	  missing_nws<-data.frame(sourceID=as.character(NA),datastream=as.character(NA)) #missing stations blank df
	  count_log_nws<-data.frame(datastream=as.character("nws"),station_count=as.numeric(0),unique=as.logical(0))
    print(paste(nws_month_filename,"MISSING empty DF made!"))
    }
print("nws pau!")

#add SCAN data
setwd(scan_daily_wd)#set data source wd
scan_month_filename<-paste0(file_date,"_scan_daily_rf.csv")
scan<-read.csv.TC(scan_month_filename)
if(!is.null(scan) && nrow(scan)>0){  #did scan month file exist? TRUE
 #subset 24hr obs
 names(scan)<-c("sourceID","date","x")
 scan$date<-as.Date(scan$date)
 scan_day<- scan[scan$date==map_date,]
 if(nrow(scan_day)>0){ #if scan_day has rows/data
  scan_day$date<-format(scan_day$date,"%Y.%m.%d")
  scan_day$x<-as.numeric(scan_day$x)
  #head(scan_day)
  scan_wide<- reshape(scan_day, idvar = "sourceID", timevar = "date", direction = "wide")
  scan_wide$datastream<-"scan"
  #head(scan_wide)
  scan_wide_merged<-merge(scan_wide,geog_meta_sub[,c("SKN","SCAN.id")],by.x="sourceID",by.y="SCAN.id",all.x=T)
  missing_scan<- scan_wide_merged[is.na(scan_wide_merged$SKN),c("sourceID","datastream")] #missing stations
  scan_wide_merged<- scan_wide_merged[!is.na( scan_wide_merged$SKN),] #remove missing stations with no SKN
  names(scan_wide_merged)[2]<-gsub("x.","X",names(scan_wide_merged)[2])#make lower case x to uppercase X for continuity
  scan_wide_merged<-scan_wide_merged[!is.na(scan_wide_merged[,2]),] #remove NA rf day vals
  #tail(scan_wide_merged)
  count_log_scan<-data.frame(datastream=as.character("scan"),station_count=as.numeric(nrow(scan_wide_merged)),unique=as.logical(0))
  print(paste(scan_month_filename,"found!",nrow(scan_wide_merged),"stations added!"))
  }else{ #else if nrow(scan_day) = 0
   scan_wide_merged<-data.frame(sourceID=as.character(NA),date_day=as.numeric(NA),datastream=as.character("scan"),SKN=as.numeric(NA))
   names(scan_wide_merged)[2]<-format(dataDate,"X%Y.%m.%d")
   missing_scan<-data.frame(sourceID=as.character(NA),datastream=as.character(NA)) #missing stations blank df
   count_log_scan<-data.frame(datastream=as.character("scan"),station_count=as.numeric(0),unique=as.logical(0))
   print(paste("NO SCAN DATA:",map_date,"!"))
   }}else{ #else if scan month df is missing make a blank df
    scan_wide_merged<-data.frame(sourceID=as.character(NA),date_day=as.numeric(NA),datastream=as.character("scan"),SKN=as.numeric(NA))
    names(scan_wide_merged)[2]<-format(dataDate,"X%Y.%m.%d")
	  missing_scan<-data.frame(sourceID=as.character(NA),datastream=as.character(NA)) #missing stations blank df
    count_log_scan<-data.frame(datastream=as.character("scan"),station_count=as.numeric(0),unique=as.logical(0))
	  print(paste(scan_month_filename,"MISSING empty DF made!"))
    }
print("scan pau!")

#add MADIS data
setwd(madis_daily_wd)#set data source wd
madis_month_filename<-paste0(file_date,"_madis_daily_rf.csv")#dynamic file name that includes month year so when month is done new file is written
madis<-read.csv.TC(madis_month_filename,HADS=TRUE)
if(!is.null(madis) && nrow(madis)>0){  #did madis month file exist? TRUE
 madis<-madis[,c("staID","date","rf")]
 names(madis)<-c("sourceID","date","x") #note 'source_id" IS "NWS.id" for madis
 madis$date<-as.Date(madis$date)
 madis<-madis[madis$data_per>=0.95,]#subset days with at least 95% data
 madis_day<-madis[madis$date==map_date,] #date sub
 if(nrow(madis_day)>0){ #if madis_day has rows/data
  madis_day$date<-format(madis_day$date,"%Y.%m.%d")
  madis_day<-madis_day[,c("sourceID","date","x")]
  madis_day$x<-as.numeric(madis_day$x)
  #head(madis_day)
  madis_wide<- reshape(madis_day, idvar = "sourceID", timevar = "date", direction = "wide")
  madis_wide$datastream<-"madis"
  #tail(madis_wide)
  madis_wide_merged_all<-merge(madis_wide,geog_meta_sub[,c("NWS.id","SKN")],by.x="sourceID",by.y="NWS.id",all.x=T)
  missing_madis<-madis_wide_merged_all[is.na(madis_wide_merged_all$SKN),c("sourceID","datastream")] #missing stations
  madis_wide_merged<-madis_wide_merged_all[!is.na(madis_wide_merged_all$SKN),] #remove missing stations with no SKN
  names(madis_wide_merged)[2]<-gsub("x.","X",names(madis_wide_merged)[2])#make lower case x to uppercase X for continuity 
  madis_wide_merged<-madis_wide_merged[!is.na(madis_wide_merged[,2]),] #remove NA rf day vals
  #tail(madis_wide_merged)
  count_log_madis<-data.frame(datastream=as.character("madis"),station_count=as.numeric(nrow(madis_wide_merged)),unique=as.logical(0)) #log of stations acquired
  print(paste(madis_month_filename,"found!",nrow(madis_wide_merged),"stations added!"))
  }else{ #else if nrow(madis_day) = 0 IE:no day data
   madis_wide_merged<-data.frame(sourceID=as.character(NA),date_day=as.numeric(NA),datastream=as.character("madis"),SKN=as.numeric(NA))
   names(madis_wide_merged)[2]<-format(dataDate,"X%Y.%m.%d")
   missing_madis<-data.frame(sourceID=as.character(NA),datastream=as.character(NA)) #missing stations blank df
   count_log_madis<-data.frame(datastream=as.character("madis"),station_count=as.numeric(0),unique=as.logical(0))#log of stations acquired
   print(paste("NO MADIS DATA:",map_date,"!"))
   }}else{ #else if madis month df is missing make a blank df
    madis_wide_merged<-data.frame(sourceID=as.character(NA),date_day=as.numeric(NA),datastream=as.character("madis"),SKN=as.numeric(NA))
    names(madis_wide_merged)[2]<-format(dataDate,"X%Y.%m.%d")
	  missing_madis<-data.frame(sourceID=as.character(NA),datastream=as.character(NA)) #missing stations blank df
	  count_log_madis<-data.frame(datastream=as.character("madis"),station_count=as.numeric(0),unique=as.logical(0))#log of stations acquired
    print(paste(madis_month_filename,"MISSING empty DF made!"))
    }
print("madis pau!")

setwd(synoptic_daily_wd)#set data source wd
synoptic_month_filename<-paste0(file_date,"_synoMeso_daily_rf.csv")#dynamic file name that includes month year so when month is done new file is written
synoptic<-read.csv.TC(synoptic_month_filename,HADS=FALSE)
if(!is.null(synoptic) && nrow(synoptic)>0){  #did synoptic month file exist? TRUE
  synoptic<-synoptic[,c("staID","date","rf")]
  names(synoptic)<-c("sourceID","date","x")
  synoptic$date<-as.Date(synoptic$date)#format as date
  synoptic<-synoptic[synoptic$data_per>=0.95,]#subset days with at least 95% data
  synoptic_day<-synoptic[synoptic$date==map_date,]
  if(nrow(synoptic_day)>0){ #if synoptic_day has rows/data
    synoptic_day$date<-format(synoptic_day$date,"%Y.%m.%d")
    synoptic_day$x<-as.numeric(synoptic_day$x)
    synoptic_wide<- reshape(synoptic_day, idvar = "sourceID", timevar = "date", direction = "wide")
    synoptic_wide$datastream<-"synoptic"
    #head(synoptic_wide)
    #tail(synoptic_wide)
    synoptic_wide_merged_all<-merge(synoptic_wide,geog_meta_sub[,c("NWS.id","SKN")],by.x="sourceID",by.y="NWS.id",all.x=T)
    missing_synoptic<-synoptic_wide_merged_all[is.na(synoptic_wide_merged_all$SKN),c("sourceID","datastream")] #missing stations
    synoptic_wide_merged<-synoptic_wide_merged_all[!is.na(synoptic_wide_merged_all$SKN),] #remove missing stations with no SKN
    names(synoptic_wide_merged)[2]<-gsub("x.","X",names(synoptic_wide_merged)[2])#make lower case x to uppercase X for continuity
    synoptic_wide_merged<-synoptic_wide_merged[!is.na(synoptic_wide_merged[,2]),] #remove NA rf day vals
    #tail(synoptic_wide_merged)
    count_log_synoptic<-data.frame(datastream=as.character("synoptic"),station_count=as.numeric(nrow(synoptic_wide_merged)),unique=as.logical(0))
    print(paste(synoptic_month_filename,"found!",nrow(synoptic_wide_merged),"stations added!"))
  }else{ #else if nrow(synoptic_day) = 0
    synoptic_wide_merged<-data.frame(sourceID=as.character(NA),date_day=as.numeric(NA),datastream=as.character("synoptic"),SKN=as.numeric(NA))
    names(synoptic_wide_merged)[2]<-format(dataDate,"X%Y.%m.%d")
    missing_synoptic<-data.frame(sourceID=as.character(NA),datastream=as.character(NA)) #missing stations blank df
    count_log_synoptic<-data.frame(datastream=as.character("synoptic"),station_count=as.numeric(0),unique=as.logical(0))
    print(paste("NO SYNOPTIC DATA:",map_date,"!"))
  }}else{ #else if synoptic month df is missing make a blank df
    synoptic_wide_merged<-data.frame(sourceID=as.character(NA),date_day=as.numeric(NA),datastream=as.character("synoptic"),SKN=as.numeric(NA))
    names(synoptic_wide_merged)[2]<-format(dataDate,"X%Y.%m.%d")
    missing_synoptic<-data.frame(sourceID=as.character(NA),datastream=as.character(NA)) #missing stations blank df
    count_log_synoptic<-data.frame(datastream=as.character("synoptic"),station_count=as.numeric(0),unique=as.logical(0))
    print(paste(synoptic_month_filename,"MISSING empty DF made!"))
  }
print("synoptic pau!")

#make and write table of all missing stations from acquired data
all_missing<-rbind(missing_hads,missing_nws,missing_scan,missing_madis,missing_synoptic)
all_missing<-all_missing[!is.na(all_missing$sourceID),] #remove no sourceID stations
if(nrow(all_missing)==0){
  all_missing<-data.frame(sourceID=NA,datastream="ALL")
}
all_missing$lastDate<-as.Date(map_date)
setwd(missing_sta_wd) #set output wd for missing station
missing_month_filename<-paste0(file_date,"_unknown_rf_sta.csv") #dynamic file name that includes month year so when month is done new file is written
print("missing file made... saving")
#conditional statement that adds obs of missing stations and removes duplicate for the month
if(file.exists(missing_month_filename)){
	rf_missing_df<-read.csv(missing_month_filename)
	rf_missing_df$lastDate<-as.Date(rf_missing_df$lastDate)
	rf_missing_df<-rbind(rf_missing_df,all_missing)
	rf_missing_df<-rf_missing_df[order(rf_missing_df$lastDate, decreasing = TRUE),]
	rf_missing_df_nodups<-rf_missing_df[!duplicated(rf_missing_df$sourceID),]
	rf_missing_df_nodups<-rf_missing_df_nodups[order(rf_missing_df_nodups$lastDate),]
	write.csv(rf_missing_df_nodups,missing_month_filename, row.names=F)
	print("monthly unknown sta table appended!")
	print(paste(missing_month_filename,"unknown sta table appended!"))
	}else{
	write.csv(all_missing,missing_month_filename, row.names=F)
	print(paste(missing_month_filename,"unknown sta table written!"))
	}
print("unknown station table below...")
print(all_missing)

#rbind all data streams and merge with geog meta
print("combinding all data...")
hads_nws_wide<-rbind.all.columns(hads_wide_merged,nws_wide_merged)
hads_nws_scan_wide<-rbind.all.columns(hads_nws_wide,scan_wide_merged)
hads_nws_scan_madis_wide<-rbind.all.columns(hads_nws_scan_wide,madis_wide_merged)
all_sta_data_wide<-rbind.all.columns(hads_nws_scan_madis_wide,synoptic_wide_merged)

print("all data combind!")

#remove stations with NA values
rf_col<-paste0("X",format(map_date,"%Y.%m.%d"))#define rf day col name
all_sta_data_wide<-all_sta_data_wide[!is.na(all_sta_data_wide[,rf_col]),] #remove na rf obs should be none

#reorder to define data stream priority
data_priority <- c("synoptic","hads","nws","madis","scan")
all_sta_data_wide<-all_sta_data_wide[order(match(all_sta_data_wide$datastream, data_priority)),] #remove dup stations by priority
dim(all_sta_data_wide)
str(all_sta_data_wide)
head(all_sta_data_wide)
tail(all_sta_data_wide)
print("combind data sorted!")

#remove dup stations based on SKN but keeping data in order defined above
all_sta_data_wide_no_dup<-all_sta_data_wide[!duplicated(all_sta_data_wide$SKN), ]
str(all_sta_data_wide_no_dup)
print(paste("station count with dups:",nrow(all_sta_data_wide)))#number of all stations
print(paste("station count without dups:",nrow(all_sta_data_wide_no_dup)))#number of unique stations

#sub cols rainfall
all_sta_data_wide_no_dup_rf<-all_sta_data_wide_no_dup[,c("SKN",rf_col)]#remove rf meta cols except SKN & RF DAY col
head(all_sta_data_wide_no_dup_rf)

#make rainfall source table and sub cols
all_sta_data_wide_no_dup_source<-all_sta_data_wide_no_dup[,c("SKN","datastream")]#remove meta cols except SKN & data stream cols
names(all_sta_data_wide_no_dup_source)<-names(all_sta_data_wide_no_dup_rf)#rename cols so source is date
head(all_sta_data_wide_no_dup_source)

#make and write table of source log station counts from acquired data
count_log_per<-rbind(count_log_hads,count_log_nws,count_log_scan)
count_log_unq<-data.frame(table(all_sta_data_wide_no_dup$datastream))
names(count_log_unq)<-c("datastream","station_count")
count_log_unq$unique<-as.logical(1)
count_log_all<-rbind(count_log_per,
				data.frame(datastream="SUBTOTAL",station_count=sum(count_log_per$station_count),unique=as.logical(0)),
				count_log_unq,
				data.frame(datastream="TOTAL",station_count=sum(count_log_unq$station_count),unique=as.logical(1)))
count_log_all$date<-map_date #add data date 
setwd(count_log_wd)
count_log_month_filename<-paste0(file_date,"_count_log_daily_rf.csv")#dynamic file name that includes month year so when month is done new file is written

#conditional statement that adds obs of per day station counts
if(file.exists(count_log_month_filename)){
	write.table(count_log_all,count_log_month_filename, row.names=F,sep = ",",col.names = F, append = T)
	print(paste(count_log_month_filename,"daily station count appended!"))
	}else{
	write.csv(count_log_all,count_log_month_filename, row.names=F)
	print(paste(count_log_month_filename,"daily station count written!"))
	}
print("final station count table below...")
print(count_log_all)

#write or append daily source data
setwd(rf_day_source_wd)#set rainfall output wd
source_month_filename<-paste0("Statewide_Daily_Source_",file_date,".csv") #dynamic file name that includes month year so when month is done new file is written

#conditional statement that adds new obs
if(file.exists(source_month_filename)){
  final_source_data <- combine_data(source_month_filename, all_sta_data_wide_no_dup_source, rf_col, geog_meta)
  write.csv(final_source_data,source_month_filename, row.names=F)
  print(paste(source_month_filename,"daily souce table appended!"))
}else{ #if month year file does not exist make a new month year file
  final_source_data<-merge(geog_meta,all_sta_data_wide_no_dup_source,by="SKN")
  write.csv(final_source_data,source_month_filename, row.names=F)
  print(paste(source_month_filename,"daily souce table written!"))
}

print("final source data table below...")
head(final_source_data)
tail(final_source_data)

#write data by creating or appending day to month for rf and source tables
#write or append daily rf data
setwd(rf_day_data_wd) #set rainfall output wd
rf_month_filename<-paste0("Statewide_Raw_Daily_RF_mm_",file_date,".csv") #dynamic file name that includes month year so when month is done new file is writen

#conditional statement that adds new obs day col
if(file.exists(rf_month_filename)){
  final_rf_data <- combine_data(rf_month_filename, all_sta_data_wide_no_dup_rf, rf_col, geog_meta)
	write.csv(final_rf_data,rf_month_filename, row.names=F)
	print(paste(rf_month_filename,"daily rainfall table appended!"))
    }else{ #if month year file does not exist make a new month year file
	final_rf_data<-merge(geog_meta,all_sta_data_wide_no_dup_rf,by="SKN")
	write.csv(final_rf_data,rf_month_filename, row.names=F)
	print(paste(rf_month_filename,"daily rainfall table written!"))
}

print("final data table below...")
head(final_rf_data)
tail(final_rf_data)

paste(dataDate,"DATA COMBIND RUN - CODE PAU!")

