#update last rf obs per station
rm(list = ls())#remove all objects in R

#packages
library(dplyr)
library(data.table)

#set dirs
mainDir <- "/home/hawaii_climate_products_container/preliminary"
codeDir<-paste0(mainDir,"/rainfall/code/source")
rf_day_data_wd <- paste0(mainDir,"/rainfall/data_outputs/tables/station_data/daily/raw/statewide") #final combine daily rainfall data
rf_last_obs_wd <- paste0(mainDir,"/rainfall/data_outputs/tables/rf_station_tracking/lastObs")

#define date
source(paste0(codeDir,"/dataDateFunc.R"))
dataDate<-dataDateMkr() #function for importing/defining date as input or as yesterday
doi<-dataDate #dataDate as date of interest (doi)
file_date<-format(doi,"%Y_%m")

#custom functions
lastObsCsvRead<-function(x){
  require(dplyr)
  require(data.table)
    df<-x
  df<-df[,c(1,grep("X",names(df)))]
  df <- as.data.frame(melt(setDT(df), id.vars = c("SKN"), variable.name = "date"))
  df <- df[!is.na(df$value),]
  df$date<-as.Date(df$date,format="X%Y.%m.%d")
  df$value<-NULL
  df <-as.data.frame(
    df %>% 
      group_by(SKN) %>%
      summarise(lastRFDate=max(date))
  )
  return(df)
}

#add master metadata with SKN and lat long
meta_url <- "https://raw.githubusercontent.com/ikewai/hawaii_wx_station_mgmt_container/main/Hawaii_Master_Station_Meta.csv"
meta<-read.csv(meta_url, colClasses=c("NESDIS.id"="character"))

#get monthly RF daily obs
setwd(rf_day_data_wd)
rf_month_filename<-paste0("Statewide_Raw_Daily_RF_mm_",file_date,".csv") #dynamic file name that includes month year so when month is done new file is writen
currentDayRF<-read.csv(rf_month_filename)

#make current month RF file obs
monthObs<-lastObsCsvRead(currentDayRF)
monthObs$lastRFDate<-as.Date(monthObs$lastRFDate)
str(monthObs)

#read last rf obs master file
setwd(rf_last_obs_wd)
lastObs<-read.csv("lastRFdayObs.csv")
lastObs$lastRFDate<-as.Date(lastObs$lastRFDate)
str(lastObs)

#combind with month last RF obs
allLastObs<-rbind(lastObs[,c("SKN","lastRFDate")],monthObs)
str(allLastObs)

#get last obs per SKN (station)
lastObsFinal<-as.data.frame(
  allLastObs %>% 
    group_by(SKN) %>%
    summarise(lastRFDate=max(lastRFDate))
)

#merge with meta and subset cols
subCols<-c("SKN","Station.Name","Observer","Network","Island","ELEV.m.","LAT","LON")
lastObsFinal<-merge(meta[subCols],lastObsFinal,by="SKN")

#overwrite old last rf obs master file
setwd(rf_last_obs_wd)
write.csv(lastObsFinal,"lastRFdayObs.csv",row.names=F)

##PAU##

