#agg daily 30yr partial gap filled rf to month-year: Aug 3 2020

rm(list = ls())#start fresh!

print("starting monthly rainfall agg")

#load packages

require(reshape2)

require(xts)



#custom functions

removeAllNA<-function(df){

	  if(length(grep("X",names(df)))>1){

		      df[rowSums(is.na(df[,grep("X",names(df))])) != ncol(df[,grep("X",names(df))]), ]

  }else{

	      df[!is.na(df[,grep("X",names(df))]),]

    }

}#remove rows where all months are NA

makeMonthlyRF<-function(rf_month_df){

	  rf_month_sub<-rf_month_df[,c(1,grep("X",names(rf_month_df)))]#keep only SKN and RF cols

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



stateSubCounty<-function(statefile,stateName,outdirCounty){

	  countList<-list(statefile$Island=="BI",statefile$Island=="MA"|statefile$Island=="KO"|statefile$Island=="MO"|statefile$Island=="LA",statefile$Island=="OA",statefile$Island=="KA")

  names(countList)<-c("BI","MN","OA","KA")

    for(j in names(countList)){

	        monCounty<-statefile[countList[[j]],]

      coDir<-paste(outdirCounty,j,sep="/")

          dir.create(coDir,showWarnings = F)

          coFileName<-paste(coDir,gsub("Statewide",j,stateName),sep="/")

	      write.csv(monCounty,coFileName,row.names = F)

	      print(paste("wrote...",coFileName))

	        }#end county loop

}#county sub function


print("seeting directories")
#dirs

inDir<-"/home/hawaii_climate_products_container/preliminary/rainfall/data_outputs/tables/station_data/daily/partial_filled/statewide"

outDir<-"/home/hawaii_climate_products_container/preliminary/rainfall/data_outputs/tables/station_data/monthly/partial_filled/statewide" #outdir of MONTHLY rf data

outdirCounty<-"/home/hawaii_climate_products_container/preliminary/rainfall/data_outputs/tables/station_data/monthly/partial_filled/county" #outdir of per county MONTHLY rf data

print(inDir)
print(outDir)
print(outdirCounty)


#add master metadata with SKN and lat long

meta_url <- "https://raw.githubusercontent.com/ikewai/hawaii_wx_station_mgmt_container/main/Hawaii_Master_Station_Meta.csv"

geo_meta<-read.csv(meta_url, colClasses=c("NESDIS.id"="character"))

head(geo_meta)



#add daily rf monthly file

setwd(inDir)

fileDate<-format(Sys.Date()-1,"%Y_%m")

fileYear<-format(Sys.Date()-1,"%Y")

print("setting file date")
print(fileDate)
print(fileYear)
setwd(inDir) #wd of & monthly and daily rf data
print("loading statewide partial filled daily rf mm file")
rf_month_df<-read.csv(paste0("Statewide_Partial_Filled_Daily_RF_mm_2022_12.csv"))

head(rf_month_df)



#make monthly

rf_month_wide<-makeMonthlyRF(rf_month_df) #subset station loop to process monthly agg 

str(rf_month_wide)

head(rf_month_wide)

tail(rf_month_wide)



#write our data file

setwd(outDir)

filename<-paste0("Statewide_Partial_Filled_Monthly_RF_mm_2022.csv")

if(file.exists(filename)){

	  yearFile<-read.csv(filename)

  yearFile<-appendMonthCol(yearDF=yearFile,monthDF=rf_month_wide,metafile=geo_meta)

    write.csv(yearFile,filename,row.names=F)

    print(paste(fileYear,"appended..."))

}else{

	write.csv(final_data_new,row.names=F)

	print(paste(fileYear,filename,"written..."))

}



#write county

stateSubCounty(yearFile,filename,outdirCounty)



#PAU
