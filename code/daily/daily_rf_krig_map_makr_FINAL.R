#this codes grabs data dates monthly file of daily rf and makes a interpolated daily rf grid map
rm(list = ls())#remove all objects in R

#set dirs
mainDir <- "/home/hawaii_climate_products_container/preliminary"
codeDir<-paste0(mainDir,"/rainfall/code/source")
meanRFwd<-"/home/hawaii_climate_products_container/preliminary/rainfall/dependencies/daily/co_daily"
varioDFwd<-"/home/hawaii_climate_products_container/preliminary/rainfall/dependencies/daily/varios"
dailyRFout<-paste0(mainDir,"/rainfall/data_outputs")

#define dates
source(paste0(codeDir,"/dataDateFunc.R"))
source(paste0(codeDir,"/daily_rf_funcs.R"))
s<-Sys.time()

dataDate<-dataDateMkr() #function for importing/defining date as input or as yesterday

##get monthly RF files partial fill & raw

# #ikewai OFF (use for dev server)
# ikeUrl<-"https://ikeauth.its.hawaii.edu/files/v2/download/public/system/ikewai-annotated-data/HCDP/workflow_data/preliminary" #url
# dailyRFname<-paste0(ikeUrl,"/rainfall/data_outputs/tables/station_data/daily/partial_filled/statewide/Statewide_Partial_Filled_Daily_RF_mm_",format(dataDate,"%Y_%m"),".csv")
# dailyRFnameRaw<-paste0(ikeUrl,"/rainfall/data_outputs/tables/station_data/daily/raw/statewide/Statewide_Raw_Daily_RF_mm_",format(dataDate,"%Y_%m"),".csv")

#local ON (use for container run) 
dailyRFname<-paste0(mainDir,"/rainfall/data_outputs/tables/station_data/daily/partial_filled/statewide/Statewide_Partial_Filled_Daily_RF_mm_",format(dataDate,"%Y_%m"),".csv")
dailyRFnameRaw<-paste0(mainDir,"/rainfall/data_outputs/tables/station_data/daily/raw/statewide/Statewide_Raw_Daily_RF_mm_",format(dataDate,"%Y_%m"),".csv")

dailyRFdf<-read.csv(dailyRFname)
head(dailyRFdf)

dailyRFdfRaw<-read.csv(dailyRFnameRaw)
head(dailyRFdfRaw)

#load variogram df data  
setwd(varioDFwd)
varioDFAll<-read.csv("all_vario_median_exp.csv")

#subset day
validOut<-dailyRFkrig(rfdailyDFmaster=dailyRFdf,
            rfdailyRawDFmaster=dailyRFdfRaw,
            varioDFAll=varioDFAll,
            data_date=dataDate,
            meanRFgridwd=meanRFwd,
            outdir=dailyRFout,
            dataVersion="experimental"
            )


print(validOut)
e<-Sys.time()
print(e-s)
print("daily rf PAU.")

# code pau
