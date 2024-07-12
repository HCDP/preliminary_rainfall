#agg daily 30yr partial gap filled rf to month-year
rm(list = ls())#start fresh!

#load packages
require(reshape2)
require(xts)

#options(error=traceback, show.error.locations = TRUE)

#set dirs
mainDir <- "/home/hawaii_climate_products_container/preliminary"
codeDir<-paste0(mainDir,"/rainfall/code/source")
inDir<-paste0(mainDir,"/rainfall/data_outputs/tables/station_data/daily/partial_filled/statewide")
outDir<-paste0(mainDir,"/rainfall/data_outputs/tables/station_data/monthly/partial_filled/statewide") #outdir of MONTHLY rf data
outDirTrack<-paste0(mainDir,"/rainfall/data_outputs/tables/rf_station_tracking/count/monthly") #outdir of MONTHLY data count
outdirCounty<-paste0(mainDir,"/rainfall/data_outputs/tables/station_data/monthly/partial_filled/county") #outdir of per county MONTHLY rf data

#define date
source(paste0(codeDir,"/dataDateFunc.R"))
dataDate<-dataDateMkr() #function for importing/defining date as input or as yesterday
fileDate<-format(dataDate,"%Y_%m")
fileYear<-format(dataDate,"%Y")
rf_col<-paste0("X",format(dataDate,"%Y.%m"))#define rf day col name
print(dataDate)

#custom functions
removeAllNA<-function(df){
  if(length(grep("X",names(df)))>1){
    df[rowSums(is.na(df[,grep("X",names(df))])) != ncol(df[,grep("X",names(df))]), ]
  }else{
    df[!is.na(df[,grep("X",names(df))]),]
  }
}#remove rows where all months are NA
makeMonthlyRF<-function(rf_day_month_df){
  dateCols<-grep("X",names(rf_day_month_df)) #date col numbers
  rf_month_sub<-rf_day_month_df[,c(1,dateCols)]#keep only SKN and RF cols
  rf_month_long<-melt(rf_month_sub, id=c("SKN"))
  stations<-unique(rf_day_month_df$SKN)
  rf_year_df<-data.frame()#blank df
  for(i in stations){
    rf_day_sta<-rf_month_long[rf_month_long$SKN==i,]
    rf_day_sta$date<-as.Date(gsub("X","",rf_day_sta$variable),format="%Y.%m.%d")
    rf_day_sta_xts<-xts(rf_day_sta$value, order.by=rf_day_sta$date)
    rf_month_sta_xts<-apply.monthly(rf_day_sta_xts, FUN=sum) #this is the magic line that aggs days to month and if a day is NA the month is NA
    rf_month_sta_df<-data.frame(SKN=i,monYr=index(rf_month_sta_xts),rf_mm=as.numeric(rf_month_sta_xts[,1]))
    rf_month_sta_df$monYr<-as.character(format(rf_month_sta_df$monYr,"X%Y.%m"))
    rf_year_df<-rbind(rf_year_df,rf_month_sta_df)
  }
  rf_month_wide<-dcast(rf_year_df, SKN ~ monYr, value.var="rf_mm")
  rf_month_wide<-removeAllNA(rf_month_wide)
  rf_month_final<-merge(geo_meta,rf_month_wide,by="SKN")
  return(rf_month_final)
}
RFMonthCheck<-function(rf_day_month_df,rf_month_df,dataDate){
  lastMonthDay<-as.Date(paste(format(dataDate,"%Y"),as.numeric(format(dataDate,"%m"))+1,01,sep="-"))-1
  allMonthDays<-seq.Date(as.Date(format(dataDate,"%Y-%m-01")),lastMonthDay,by="days")
  nMonDays<-length(allMonthDays)
  countyMonthDayRF<-list(rf_day_month_df$Island=="BI",rf_day_month_df$Island=="MA"|rf_day_month_df$Island=="KO"|rf_day_month_df$Island=="MO"|rf_day_month_df$Island=="LA",rf_day_month_df$Island=="OA",rf_day_month_df$Island=="KA")
  counties<-c("BI","MN","OA","KA")
  names(countyMonthDayRF)<-c("BI","MN","OA","KA")
  StaCounties<-ifelse(rf_month_df$Island=="MA"|rf_month_df$Island=="KO"|rf_month_df$Island=="MO"|rf_month_df$Island=="LA","MN",rf_month_df$Island)
  allDataCheck<-data.frame()
  for(c in counties){
    co_rf_day_month_df<-rf_day_month_df[countyMonthDayRF[[c]],]
    co_rf_month_df<-rf_month_df[StaCounties==c,]
    dateCols<-grep("X",names(co_rf_day_month_df)) #date col numbers
    #co_rf_day_month_df<-co_rf_day_month_df[,dateCols]
    #co_rf_day_month_df <- co_rf_day_month_df[,colSums(is.na(co_rf_day_month_df))<nrow(co_rf_day_month_df)]    
    colDates<-as.Date(names(co_rf_day_month_df)[dateCols],format="X%Y.%m.%d")
    co_rf_day_month_df_dates<-co_rf_day_month_df[,c(1,dateCols)]
    ccRFdates<-sum(complete.cases(co_rf_day_month_df_dates))
    CompleteMonthRFstations<-ifelse(length(grep("X",names(co_rf_day_month_df_dates)))==nMonDays,ccRFdates,0)
    MonthyPCstations<-nrow(co_rf_month_df)-ccRFdates
    coDataCheck<-data.frame(yearMonth=format(dataDate,"X%Y.%m"),County=c,
                            MaximumRFstation=nrow(co_rf_day_month_df),
                            CompleteMonthRFstations=CompleteMonthRFstations,
                            MonthyPCstations=MonthyPCstations,
                            TotalMonthRFstations=MonthyPCstations+CompleteMonthRFstations,
                            DailyQAQCperformed=MonthyPCstations==0,
                            TotalMissingRFdays=sum(!allMonthDays%in%colDates),  
                            MissingRFdays=paste(allMonthDays[!allMonthDays%in%colDates],collapse = ","))
    allDataCheck<-rbind(allDataCheck,coDataCheck)
  }
  dateColsAll<-grep("X",names(rf_day_month_df)) #date col numbers
  colDatesAll<-as.Date(names(rf_day_month_df)[dateCols],format="X%Y.%m.%d")
  allDataCheck<-rbind(allDataCheck,
                      data.frame(yearMonth=format(dataDate,"X%Y.%m"),County="All",
                                 MaximumRFstation=sum(allDataCheck$MaximumRFstation),
                                 CompleteMonthRFstations=sum(allDataCheck$CompleteMonthRFstations),
                                 MonthyPCstations=sum(allDataCheck$MonthyPCstations),
                                 TotalMonthRFstations=sum(allDataCheck$TotalMonthRFstations),
                                 DailyQAQCperformed=as.logical(min(allDataCheck$DailyQAQCperformed)),
                                 TotalMissingRFdays=sum(!allMonthDays%in%colDatesAll),    
                                 MissingRFdays=paste(allMonthDays[!allMonthDays%in%colDatesAll],collapse = ","))
                      )
  
  return(allDataCheck)
}
MonthPC2rf<-function(doi,missingCo){
  
  firstDate<-as.Date(format(doi,"%Y-%m-01"))
  lastDate<-doi+1
  
  dtoi1<-strptime(paste(firstDate, "00:00:00"),format="%Y-%m-%d %H:%M:%S")
  dtoi2<-strptime(paste(lastDate, "00:00:00"),format="%Y-%m-%d %H:%M:%S")
  attr(dtoi1,"tzone") <- "Pacific/Honolulu" #convert TZ attribute to HST
  attr(dtoi2,"tzone") <- "Pacific/Honolulu" #convert TZ attribute to HST
  
  #hads url
  URL<-"https://ikeauth.its.hawaii.edu/files/v2/download/public/system/ikewai-annotated-data/HCDP/workflow_data/preliminary/data_aqs/data_outputs/hads/parse/"
  
  #get first day of month hads
  file<-paste0(format(firstDate,"%Y%m%d"),"_hads_parsed.csv")
  firstHADS<-read.csv(paste0(URL,file))
  
  #get last day of month hads
  file<-paste0(format(lastDate,"%Y%m%d"),"_hads_parsed.csv")
  lastHADS<-read.csv(paste0(URL,file))
  
  #rbind together
  allHADS<-rbind(firstHADS,lastHADS)
  hads_pc<-subset(allHADS,var=="PC")# subset precip only
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
  names(hads_pc1)[ncol(hads_pc1)]<-as.character(dtoi1)
  row.names(hads_pc1)<-NULL
  #pc2
  hads_pc2<-hads_pc[hads_pc$obs_time==dtoi2,]  #subset only matching date time pc obs
  hads_pc2<-hads_pc2[,c("staID","value")]
  names(hads_pc2)[ncol(hads_pc2)]<-as.character(dtoi2)
  row.names(hads_pc2)<-NULL
  
  #merge
  mergedPC<-merge(hads_pc1,hads_pc2,by="staID")
  mergedPC$PCdiff<-mergedPC[,as.character(dtoi2)]-mergedPC[,as.character(dtoi1)]
  mergedPC<-mergedPC[mergedPC$PCdiff>=0 & !is.na(mergedPC$PCdiff),c("staID","PCdiff")] #subset minus values
  
  #add metadata
  meta_mergedPC<-merge(geo_meta,mergedPC,by.x="NESDIS.id",by.y="staID")
  subCols<-c(names(geo_meta),"PCdiff")
  
  meta_mergedPC$StaCounties<-ifelse(meta_mergedPC$Island=="MA"|meta_mergedPC$Island=="KO"|meta_mergedPC$Island=="MO"|meta_mergedPC$Island=="LA","MN",
                      meta_mergedPC$Island)
  #make month matching DF 
  meta_mergedPC<-meta_mergedPC[meta_mergedPC$StaCounties %in% missingCo,subCols]
  names(meta_mergedPC)[ncol(meta_mergedPC)]<-format(firstDate,"X%Y.%m") #change RF name
  row.names(meta_mergedPC)<-NULL
  return(meta_mergedPC)
}#end MonthPC2rf
appendMonthCol<-function(yearDF,monthDF,metafile,rf_col){
  yearDFcols<-names(yearDF)
  yearDFcols<-yearDFcols[yearDFcols!=rf_col]
  sub_cols<-c(1,grep("X",yearDFcols))
  if(length(sub_cols)>1){
    yearDFsub<-yearDF[,sub_cols]#keep only SKN and monthly RF cols
    monthDF<-monthDF[,c(1,grep("X",names(monthDF)))]#keep only SKN and monthly RF cols
    yearDFsub<-merge(metafile,yearDFsub,by="SKN",all=T)
    yearDFsub<-yearDFsub[,c(1,grep("X",names(yearDFsub)))]#keep only SKN and monthly RF cols
    monthDF<-merge(metafile,monthDF,by="SKN",all=T)
    monthDF<-monthDF[,c(1,grep("X",names(monthDF)))]#keep only SKN and monthly RF cols
    yearDFsub<-merge(yearDFsub,monthDF,by="SKN")
    yearDFsub<-removeAllNA(yearDFsub)
    yearFinal<-merge(metafile,yearDFsub,by="SKN")
  }else{
    yearFinal<-monthDF
  }
  message("month added to year!")
  return(yearFinal)
}
stateSubCounty<-function(stateFile,stateName=NA,outdirCounty=NA,writeCo=F){
  countList<-list(stateFile$Island=="BI",stateFile$Island=="MA"|stateFile$Island=="KO"|stateFile$Island=="MO"|stateFile$Island=="LA",stateFile$Island=="OA",stateFile$Island=="KA")
  names(countList)<-c("BI","MN","OA","KA")
  stateCoList<-list()
  for(j in names(countList)){
    monCounty<-stateFile[countList[[j]],]
    stateCoList[[j]]<-monCounty
    if(writeCo){
      coDir<-paste(outdirCounty,j,sep="/")
      dir.create(coDir,showWarnings = F)
      coFileName<-paste(coDir,gsub("Statewide",j,stateName),sep="/")
      write.csv(monCounty,coFileName,row.names = F)
      message(paste("wrote...",coFileName))
    }
  }#end county loop
  return(stateCoList)
}#county sub function
rbind.all.columns <- function(x, y) {     #function to smart rbind
  x.diff <- setdiff(colnames(x), colnames(y))
  y.diff <- setdiff(colnames(y), colnames(x))
  x[, c(as.character(y.diff))] <- NA 
  y[, c(as.character(x.diff))] <- NA 
  return(rbind(x, y))}
#add master metadata with SKN and lat long
meta_url <- "https://raw.githubusercontent.com/ikewai/hawaii_wx_station_mgmt_container/main/Hawaii_Master_Station_Meta.csv"
geo_meta<-read.csv(meta_url, colClasses=c("NESDIS.id"="character"))
#head(geo_meta)

#add daily rf monthly file local
setwd(inDir) #wd of & monthly and daily rf data
rf_month_df<-read.csv(paste0("Statewide_Partial_Filled_Daily_RF_mm_",fileDate,".csv"))
head(rf_month_df)

# #add daily rf monthly file url
# URL<-"https://ikeauth.its.hawaii.edu/files/v2/download/public/system/ikewai-annotated-data/HCDP/workflow_data/preliminary/rainfall/data_outputs/tables/station_data/daily/partial_filled/statewide/"
# rf_month_df<-read.csv(paste0(URL,"Statewide_Partial_Filled_Daily_RF_mm_",fileDate,".csv"))
# head(rf_month_df)

#make monthly
rf_month_wide<-makeMonthlyRF(rf_day_month_df=rf_month_df) #subset station loop to process monthly agg 
#str(rf_month_wide)
#head(rf_month_wide)
#tail(rf_month_wide)

#monthly rf check 
rf_month_track<-RFMonthCheck(rf_day_month_df=rf_month_df,
                             rf_month_df=rf_month_wide,
                             dataDate=dataDate)
print("station tracking...")
print(rf_month_track)

#check each county has station data and run stop gap month hads process if county has no stations
nSta<-10 #min stations per county
missingCo<-rf_month_track[rf_month_track$CompleteMonthRFstations<nSta & rf_month_track$County!="All","County"]

#get hads PC stations for county month values
if(length(missingCo)>0){
  missingCoRF<-MonthPC2rf(doi=dataDate,missingCo=missingCo)
  rf_month_wide<-rbind(rf_month_wide,missingCoRF)
  rf_month_wide<-rf_month_wide[!duplicated(rf_month_wide$SKN),]
  rf_month_track<-RFMonthCheck(rf_month_df,rf_month_wide,dataDate)
  print("station tracking REDO...")
  print(rf_month_track)
}

#check data availible to krig each county
coStaCheck<-rf_month_track[rf_month_track$County!="All","TotalMonthRFstations"]>3
lowCounty<-paste(rf_month_track[coStaCheck,"County"],collapse=" &")

if(!all(coStaCheck)){
  stop(paste(lowCounty,"did not have enough stations to perform monthly RF krigging!!"))
}else{ #write files and finish script
    
  #save monthly rf annual file table
  setwd(outDir)
  filename<-paste0("Statewide_Partial_Filled_Monthly_RF_mm_",fileYear,".csv") #dynamic file name that includes year of file

  #append or write new annual monthly rf file
  filename<-paste0("Statewide_Partial_Filled_Monthly_RF_mm_",fileYear,".csv")
  if(file.exists(filename)){ #check if downloaded file is in wd
    yearFile<-read.csv(filename)
    yearFile<-appendMonthCol(yearDF=yearFile,monthDF=rf_month_wide,metafile=geo_meta,rf_col=rf_col)
    write.csv(yearFile,filename,row.names=F)
    print(paste(fileYear,"appended..."))
  }else{ #if file did not exist/download write new file
    yearFile<-rf_month_wide
    write.csv(yearFile,filename,row.names=F)
    print(paste(fileYear,filename,"written..."))
  }
  
  #sub county data and save
  stateCoList<-stateSubCounty(stateFile=yearFile,stateName=filename,outdirCounty=outdirCounty,writeCo=TRUE)
  #head(stateCoList) #county sub list
  
  #write monthly check
  setwd(outDirTrack)
  rf_month_track_filename<-paste0(fileYear,"_count_log_monthly_rf.csv") #dynamic file name
  
  if(file.exists(rf_month_track_filename)){
    rf_month_track<-rbind.all.columns(read.csv(rf_month_track_filename),rf_month_track) #load old table and append new table
    write.csv(rf_month_track,rf_month_track_filename, row.names=F) #write appended table
    print(paste(rf_month_track_filename,"monthly station count appended!"))
  }else{
    write.csv(rf_month_track,rf_month_track_filename, row.names=F)
    print(paste(rf_month_track_filename,"monthly station count written!"))
  }
  print("day to month rf pau!")
}
#PAU