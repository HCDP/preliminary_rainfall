#agg daily 30yr partial gap filled rf to month-year
rm(list = ls())#start fresh!

#load packages
require(reshape2)
require(xts)

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
print(dataDate)

#custom functions
removeAllNA<-function(df){
  if(length(grep("X",names(df)))>1){
    df[rowSums(is.na(df[,grep("X",names(df))])) != ncol(df[,grep("X",names(df))]), ]
  }else{
    df[!is.na(df[,grep("X",names(df))]),]
  }
}#remove rows where all months are NA
makeMonthlyRF<-function(rf_month_df){
  dateCols<-grep("X",names(rf_month_df)) #date col numbers
  rf_month_sub<-rf_month_df[,c(1,dateCols)]#keep only SKN and RF cols
  rf_month_long<-melt(rf_month_sub, id=c("SKN"))
  stations<-unique(rf_month_df$SKN)
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
RFMonthCheck<-function(rf_month_df,dataDate){
  allMonthDays<-seq.Date(as.Date(format(dataDate,"%Y-%m-01")),dataDate,by="days")
  countyMonthDayRF<-list(rf_month_df$Island=="BI",rf_month_df$Island=="MA"|rf_month_df$Island=="KO"|rf_month_df$Island=="MO"|rf_month_df$Island=="LA",rf_month_df$Island=="OA",rf_month_df$Island=="KA")
  names(countyMonthDayRF)<-c("BI","MN","OA","KA")
  allDataCheck<-data.frame()
  counties<-c("BI","MN","OA","KA")
  for(c in counties){
    co_rf_month_df<-rf_month_df[countyMonthDayRF[[c]],]
    dateCols<-grep("X",names(co_rf_month_df)) #date col numbers
    #co_rf_month_df<-co_rf_month_df[,dateCols]
    #co_rf_month_df <- co_rf_month_df[,colSums(is.na(co_rf_month_df))<nrow(co_rf_month_df)]    
    colDates<-as.Date(names(co_rf_month_df)[dateCols],format="X%Y.%m.%d")
    coDataCheck<-data.frame(yearMonth=format(dataDate,"X%Y.%m"),County=c,
                            CountRFstation=nrow(co_rf_month_df),
                            TotalMissingRFdays=sum(!allMonthDays%in%colDates),    
                            MissingRFdays=paste(allMonthDays[!allMonthDays%in%colDates],collapse = ","))
    allDataCheck<-rbind(allDataCheck,coDataCheck)
  }
  dateColsAll<-grep("X",names(rf_month_df)) #date col numbers
  colDatesAll<-as.Date(names(rf_month_df)[dateCols],format="X%Y.%m.%d")
  allDataCheck<-rbind(allDataCheck,
                      data.frame(yearMonth=format(dataDate,"X%Y.%m"),County="All",
                                 CountRFstation=sum(allDataCheck$CountRFstation),
                                 TotalMissingRFdays=sum(!allMonthDays%in%colDatesAll),    
                                 MissingRFdays=paste(allMonthDays[!allMonthDays%in%colDatesAll],collapse = ","))
                      )
  return(allDataCheck)
}
appendMonthCol<-function(yearDF,monthDF,metafile){
  yearDFsub<-yearDF[,c(1,grep("X",names(yearDF)))]#keep only SKN and monthly RF cols
  monthDF<-monthDF[,c(1,grep("X",names(monthDF)))]#keep only SKN and monthly RF cols
  yearDFsub<-merge(metafile,yearDFsub,by="SKN",all=T)
  yearDFsub<-yearDFsub[,c(1,grep("X",names(yearDFsub)))]#keep only SKN and monthly RF cols
  monthDF<-merge(metafile,monthDF,by="SKN",all=T)
  monthDF<-monthDF[,c(1,grep("X",names(monthDF)))]#keep only SKN and monthly RF cols
  if(sum(names(yearDFsub) %in% names(monthDF))>1){
    message("dup month! NOT added")
    return(yearDF)
  }else{
    yearDFsub<-merge(yearDFsub,monthDF,by="SKN")
    yearDFsub<-removeAllNA(yearDFsub)
    yearFinal<-merge(metafile,yearDFsub,by="SKN")
    message("month added to year!")
    return(yearFinal)
    }
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

#add master metadata with SKN and lat long
meta_url <- "https://raw.githubusercontent.com/ikewai/hawaii_wx_station_mgmt_container/main/Hawaii_Master_Station_Meta.csv"
geo_meta<-read.csv(meta_url, colClasses=c("NESDIS.id"="character"))
head(geo_meta)

#add daily rf monthly file
setwd(inDir) #wd of & monthly and daily rf data
rf_month_df<-read.csv(paste0("Statewide_Partial_Filled_Daily_RF_mm_",fileDate,".csv"))
head(rf_month_df)

#make monthly
rf_month_wide<-makeMonthlyRF(rf_month_df) #subset station loop to process monthly agg 
str(rf_month_wide)
head(rf_month_wide)
tail(rf_month_wide)

#save monthly rf annual file table
setwd(outDir)
filename<-paste0("Statewide_Partial_Filled_Monthly_RF_mm_",fileYear,".csv") #dynamic file name that includes year of file

#append or write new annual monthly rf file
filename<-paste0("Statewide_Partial_Filled_Monthly_RF_mm_",fileYear,".csv")
if(file.exists(filename)){ #check if downloaded file is in wd
  yearFile<-read.csv(filename)
  yearFile<-appendMonthCol(yearDF=yearFile,monthDF=rf_month_wide,metafile=geo_meta)
  write.csv(yearFile,filename,row.names=F)
  print(paste(fileYear,"appended..."))
}else{ #if file did not exist/download write new file
  yearFile<-rf_month_wide
  write.csv(yearFile,filename,row.names=F)
  print(paste(fileYear,filename,"written..."))
}

#sub county data and save
stateCoList<-stateSubCounty(stateFile=yearFile,stateName=filename,outdirCounty=outdirCounty,writeCo=TRUE)
head(stateCoList) #county sub list

#monthly rf check 
rf_month_track<-RFMonthCheck(rf_month_df,dataDate)

#write monthly check
setwd(outDirTrack)
rf_month_track_filename<-paste0(fileYear,"_count_log_monthly_rf.csv") #dynamic file name

if(file.exists(rf_month_track_filename)){
  write.table(rf_month_track,rf_month_track_filename, row.names=F,sep = ",",col.names = F, append = T)
  print(paste(rf_month_track_filename,"monthly station count appended!"))
}else{
  write.csv(rf_month_track,rf_month_track_filename, row.names=F)
  print(paste(rf_month_track_filename,"monthly station count written!"))
}
print("final station count table below...")
print(rf_month_track) #station count and missing day check per county and statewide

#PAU
