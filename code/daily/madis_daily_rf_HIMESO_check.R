#this code grabs all raw madis data for HI from gateway
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
dataDate<-dataDateMkr("2023-05-25") #function for importing/defining date as input or as yesterday
currentDate<-dataDate #dataDate as currentDate (yesterday)

#read MADIS parsed table from ikewai data portal
ikeUrl<-"https://ikeauth.its.hawaii.edu/files/v2/download/public/system/ikewai-annotated-data/HCDP/workflow_data/preliminary_test" #url
madis_filename<-paste0(format((currentDate),"%Y%m%d"),"_madis_parsed.csv") #dynamic file name that includes date
all_madis<-read.csv(paste0(ikeUrl,"/data_aqs/data_outputs/madis/parse/",madis_filename))
head(all_madis)
tail(all_madis)
unique(all_madis$varname)

hiMesonet<-c('017HI','016HI','001HI','019HI','002HI','013HI','003HI','020HI','023HI','022HI','021HI','018HI','029HI','025HI','027HI','005HI','006HI','007HI','008HI','009HI','010HI','011HI','012HI','014HI','015HI','024HI')
#hiMesonet<-c('021HI','022HI','027HI','025HI','008HI','018HI','002HI','020HI','017HI','001HI','016HI','028HI','019HI','011HI','012HI','026HI','015HI','014HI','024HI')
hiMeso<-all_madis[all_madis$stationId %in% hiMesonet,]
unique(hiMeso$varname)
