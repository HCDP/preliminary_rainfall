#this codes grabs data dates montlth file of daily rf and makes a interpolated daily grid
rm(list = ls())#remove all objects in R

#set dirs
mainDir <- "/home/hawaii_climate_products_container/preliminary"
codeDir<-paste0(mainDir,"/rainfall/code/source")
meanRFwd<-"/home/hawaii_climate_products_container/preliminary/rainfall/dependencies/daily/co_daily"
varioDFwd<-"/home/hawaii_climate_products_container/preliminary/rainfall/dependencies/daily/varios"
dailyRFout<-paste0(mainDir,"/rainfall/dailyRF_test")

#define dates
source(paste0(codeDir,"/dataDateFunc.R"))
source(paste0(codeDir,"/daily_rf_funcs.R"))
s<-Sys.time()

dataDate<-dataDateMkr() #function for importing/defining date as input or as yesterday

##get monthly RF file
#partial fill 
ikeUrl<-"https://ikeauth.its.hawaii.edu/files/v2/download/public/system/ikewai-annotated-data/HCDP/workflow_data/preliminary" #url
dailyRFname<-paste0(ikeUrl,"/rainfall/data_outputs/tables/station_data/daily/partial_filled/statewide/Statewide_Partial_Filled_Daily_RF_mm_",format(dataDate,"%Y_%m"),".csv")
dailyRFdf<-read.csv(dailyRFname)
head(dailyRFdf)

#raw
dailyRFnameRaw<-paste0(ikeUrl,"/rainfall/data_outputs/tables/station_data/daily/raw/statewide/Statewide_Raw_Daily_RF_mm_",format(dataDate,"%Y_%m"),".csv")
dailyRFdfRaw<-read.csv(dailyRFnameRaw)
head(dailyRFdfRaw)

#load variogram df data  
setwd(varioDFwd)
varioDFAll<-read.csv("all_vario_median_exp.csv")

#subset day
testOut<-dailyRFkrig(rfdailyDFmaster=dailyRFdf,
            rfdailyRawDFmaster=dailyRFdfRaw,
            varioDFAll=varioDFAll,
            data_date=dataDate,
            meanRFgridwd=meanRFwd,
            outdir=dailyRFout,
            dataVersion="preliminary"
            )


print(testOut)
e<-Sys.time()
e-s
#pau
