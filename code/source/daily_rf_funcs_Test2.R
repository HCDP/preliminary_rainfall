library(automap)
library(raster)
library(gstat)
library(hydroGOF)
library(Metrics)
library(doParallel)
library(foreach)
library(lubridate)
library(parallel)


#custom funcs
#functions
rbindAll <- function(x, y) {     #function to smart rbind
  x.diff <- setdiff(colnames(x), colnames(y))
  y.diff <- setdiff(colnames(y), colnames(x))
  x[, c(as.character(y.diff))] <- NA 
  y[, c(as.character(x.diff))] <- NA 
  return(rbind(x, y))}

scale_values <- function(x, ...) {(x - min(x, ...)) / (max(x, ...) - min(x, ...))}

rankMean<-function(x){
  scaleDF<-data.frame(
    RMSE=scale_values(x$rmse_rf_mm),
    MAE=scale_values(x$mae_rf_mm),
    BIAS=scale_values(abs(x$bias_rf_mm)))
  meanScore<-rowMeans(scaleDF)
  return(meanScore)
}

realFill<-function(all,real){return(all$SKN %in% real$SKN)}

getValidMets<-function(county=NA,data_date,addC=NA,stationCount,useVario=NA,nugFixZero=NA,vario=NA,RF_day_raw=NA,loocv_df,statewide=F){
  if(statewide){
    if(max(loocv_df$obs_rf)==0){
      validation<-data.frame(
        statewide=TRUE,
        counties=paste(c("BI","MN","OA","KA"),collapse = " "),
        date=data_date,
        stationCount=stationCount,
        stationCountReal=length(unique(loocv_df$SKN)),
        stationCountFill=stationCount-length(unique(loocv_df$SKN)),
        staRFmmMin=min(loocv_df$obs_rf),
        staRFmmMax=max(loocv_df$obs_rf),
        noRF=TRUE,
        rsq_rf_mm=NA,rmse_rf_mm=NA,mae_rf_mm=NA,
        bias_rf_mm=NA,nse_rf_mm=NA,kge_rf_mm=NA,
        rsq_rf_anom=NA,rmse_rf_anom=NA,mae_rf_anom=NA,
        bias_rf_anom=NA,nse_rf_anom=NA,kge_rf_anom=NA)
    }else{
      validation<-data.frame(
        statewide=TRUE,
        counties=paste(c("BI","MN","OA","KA"),collapse = " "),
        date=data_date,
        stationCount=stationCount,
        stationCountReal=length(unique(loocv_df$SKN)),
        stationCountFill=stationCount-length(unique(loocv_df$SKN)),
        staRFmmMin=min(loocv_df$obs_rf),
        staRFmmMax=max(loocv_df$obs_rf),
        noRF=FALSE,
        rsq_rf_mm=summary(lm(loocv_df$obs_rf~loocv_df$pred_rf))$r.squared,
        rmse_rf_mm=Metrics::rmse(loocv_df$obs_rf,loocv_df$pred_rf),
        mae_rf_mm=Metrics::mae(loocv_df$obs_rf,loocv_df$pred_rf),
        bias_rf_mm=bias(loocv_df$obs_rf,loocv_df$pred_rf),
        nse_rf_mm=NSE(sim=loocv_df$pred_rf,obs=loocv_df$obs_rf),
        kge_rf_mm=KGE(sim=loocv_df$pred_rf,obs=loocv_df$obs_rf,out.type="single"),
        rsq_rf_anom=summary(lm(loocv_df$obs_anom~loocv_df$pred_anom))$r.squared,
        rmse_rf_anom=Metrics::rmse(loocv_df$obs_anom,loocv_df$pred_anom),
        mae_rf_anom=Metrics::mae(loocv_df$obs_anom,loocv_df$pred_anom),
        bias_rf_anom=bias(loocv_df$obs_anom,loocv_df$pred_anom), 
        nse_rf_anom=NSE(sim=loocv_df$pred_anom,obs=loocv_df$obs_anom),
        kge_rf_anom=KGE(sim=loocv_df$pred_anom,obs=loocv_df$obs_anom,out.type="single")
      )}#rainfall recorded
  }else{ #per county
    validation<-data.frame(
      county=county,
      date=data_date,
      stationCount=stationCount,
      stationCountReal=nrow(RF_day_raw),
      stationCountFill=stationCount-nrow(RF_day_raw),
      staRFmmMin=min(loocv_df$obs_rf),
      staRFmmMax=max(loocv_df$obs_rf),
      noRF=FALSE,
      addC=addC,
      fixedVario=useVario,
      nugFixZero=nugFixZero,
      mod=as.character(vario$var_model$model[2]),
      nugget=vario$var_model[1,2],
      range=vario$var_model[2,3],
      sill=vario$var_model[2,2],
      rsq_rf_mm=summary(lm(loocv_df$obs_rf~loocv_df$pred_rf))$r.squared,
      rmse_rf_mm=Metrics::rmse(loocv_df$obs_rf,loocv_df$pred_rf),
      mae_rf_mm=Metrics::mae(loocv_df$obs_rf,loocv_df$pred_rf),
      bias_rf_mm=bias(loocv_df$obs_rf,loocv_df$pred_rf),
      nse_rf_mm=NSE(sim=loocv_df$pred_rf,obs=loocv_df$obs_rf),
      kge_rf_mm=KGE(sim=loocv_df$pred_rf,obs=loocv_df$obs_rf,out.type="single"),
      rsq_rf_anom=summary(lm(loocv_df$obs_anom~loocv_df$pred_anom))$r.squared,
      rmse_rf_anom=Metrics::rmse(loocv_df$obs_anom,loocv_df$pred_anom),
      mae_rf_anom=Metrics::mae(loocv_df$obs_anom,loocv_df$pred_anom),
      bias_rf_anom=bias(loocv_df$obs_anom,loocv_df$pred_anom), 
      nse_rf_anom=NSE(sim=loocv_df$pred_anom,obs=loocv_df$obs_anom),
      kge_rf_anom=KGE(sim=loocv_df$pred_anom,obs=loocv_df$obs_anom,out.type="single")
    )}
  return(validation)
} #end validation metrics func

metamaker<-function(map_validation_df,grid,filenames,datatype,rfDay,statewide=F,state_validation=NA,loocv_df=NA){
  validRanks<-c(0.025,0.05,0.075,0.1) #ranks of high,good,moderate,low respectively 
  map_validation_df$qualMetric<-map_validation_df$mae_rf_mm/map_validation_df$staRFmmMax
  
  #packages
  require(raster)
  #get dates
  dataStartDate<-rfDay
  dataEndDate<-rfDay+1
  dataconvert<-function(x){return(ifelse(is.numeric(x),as.character(round(x,5)),as.character(x)))}
  function(x){return(ifelse(is.numeric(x),as.character(round(x,5)),as.character(x)))}
  if(statewide==F){
    valid_meta<-data.frame(attribute=as.character(names(map_validation_df)),value=as.character(lapply(map_validation_df[1,], dataconvert)))
    
    countyText<-ifelse(map_validation_df$county=="ka"|map_validation_df$county=="KA","Kauai County",
                       ifelse(map_validation_df$county=="oa"|map_validation_df$county=="OA","Honolulu County (Oahu)",
                              ifelse(map_validation_df$county=="mn"|map_validation_df$county=="MN","Maui County (Maui, Lanai, Molokai & Kahoolawe)",
                                     if(map_validation_df$county=="bi"|map_validation_df$county=="BI"){"Hawaii county"})))
    
    co_statement<-as.character(paste("This",format(dataStartDate,"%B %d %Y"),"rainfall map of",countyText, 
                                     "is a high spatial resolution (~250m) gridded prediction of cumulative rainfall in millimeters from midnight HST",
                                     format(dataStartDate,"%b %d %Y"),"to midnight HST",format(dataEndDate,"%b %d %Y."),
                                     ifelse(map_validation_df$noRF,as.character(paste("All",map_validation_df$stationCount, "observed and statistically gap filled unique station locations within",countyText, "county reported no measured or estimated rainfall. Without any observed rainfall the county wide map was created by setting all pixels to 0. As such, no validation metric data is availible, and an rainfall SE map was not produced.", 
                                                                                      "Unknown rainfall errors could be present, if unreported rainfall occured at locations without known station observations. If so this map might be underestimating accual rainfall in unknown areas. All maps are subject to change as new data becomes available or unknown errors are corrected in reoccurring versions.")),
                                            #IDW statement here
                                            as.character(paste("This was produced using a climate-aided modified automatic kriging interpolation of a log transformed daily rainfall anomaly ratio calculated from a mean monthly daily disaggrigated rf map, plus a consteint (observed mm + 1 / (mean daily mm + c). This kriging process used", 
                                                               map_validation_df$stationCount, "unique station locations within",countyText,
                                                               "and their",format(dataStartDate,"%B %Y"), 
                                                               "recorded and/or estimated rainfall (mm) totals. A leave one out cross validation (LOOCV) between observed station data and the interpolated estimate used in this map produced a mean absolute error (MAE) of:",
                                                               round(map_validation_df$mae_rf_mm,3),"with a maximum observed rainfall of",map_validation_df$staRFmmMax,"mm, meaning this",format(dataStartDate,"%b %d %Y"),
                                                               countyText,"daily rainfall (mm) map is a", 
                                                               ifelse(map_validation_df$qualMetric>=validRanks[1],"high quality estimate of daily rainfall.",
                                                                      ifelse(map_validation_df$qualMetric>=validRanks[2],"good quality estimate of daily rainfall.",
                                                                             ifelse(map_validation_df$qualMetric>=validRanks[3],"moderate quality estimate of daily rainfall.",
                                                                                    ifelse(map_validation_df$qualMetric>=validRanks[4],"low quality estimate of daily rainfall, and should be used with dilligence.","lowest quality estimate of daily rainfall, and should be used with caution.")))),
                                                               "All maps are subject to change as new data becomes available or unknown errors are corrected in reoccurring versions.", 
                                                               "Errors in rainfall estimates do vary over space meaning any gridded rainfall value, even on higher quality maps, could still produce incorrect estimates.",
                                                               "Check standard error (SE) maps to better understand spatial estimates of prediction error")))))
    
    keywords<-as.character(paste(ifelse(map_validation_df$county=="ka"|map_validation_df$county=="KA","Kauai,",
                                        ifelse(map_validation_df$county=="oa"|map_validation_df$county=="OA","Oahu,",
                                               ifelse(map_validation_df$county=="mn"|map_validation_df$county=="MN","Maui, Lanai, Molokai, Kahoolawe,)",
                                                      if(map_validation_df$county=="bi"|map_validation_df$county=="BI"){"Hawaii Island,"}))),"Hawaiian Islands, rainfall prediction, daily precipitation, rainfall, climate, spatial interpolation, kriging"))
    
    file_meta<-data.frame(dataStatement=co_statement,
                          keywords=keywords,
                          county=map_validation_df$county,
                          dataStartDate=dataStartDate,
                          dataEndDate=dataEndDate,
                          dateProduced=Sys.Date(),
                          dataVersionType=datatype,
                          RFstationFile=filenames[1],
                          RFmmGridFile=filenames[2],
                          RFmmSEgridFile=filenames[3],
                          GeoCoordUnits="Decimal Degrees",
                          GeoCoordRefSystem=as.character(crs(grid)),
                          Xresolution=xres(grid),
                          Yresolution=yres(grid),
                          ExtentXmin=extent(grid)[1],
                          ExtentXmax=extent(grid)[2],
                          ExtentYmin=extent(grid)[3],
                          ExtentYmax=extent(grid)[4]
    )
    
    final_meta<-rbind(data.frame(attribute=as.character(names(file_meta)),value=as.character(lapply(file_meta[1,], dataconvert))),
                      valid_meta[-c(1,2),],
                      data.frame(attribute=as.character(c("credits", "contacts")),
                                 value=as.character(c("All data produced by University of Hawaii at Manoa Water Resource Research Center (WRRC). Support provided by Change Hawaii EPSCoR funded by the National Science Foundation under EPSCoR Research Infrastructure Improvement Award #OIA-2149133",
                                                      "Matthew Lucas (mplucas@hawaii.edu), Keri Kodama (kodamak8@hawaii.edu), Ryan Longman (rlongman@hawaii.edu), Thomas Giambelluca (thomas@hawaii.edu)")))		
    )
    final_meta$attribute<-as.character(final_meta$attribute)
    row.names(final_meta)<-NULL
    return(final_meta)
  }# county meta
  
  if(statewide){
    #make state per county validation metrics
    map_validation_df_t<-as.data.frame(t(map_validation_df))
    colapseValidation_t <- as.character(apply( map_validation_df_t , 1 , paste , collapse = ", " ))
    colapseValidation_t<-colapseValidation_t[-2]
    rf_validation<-as.data.frame(t(colapseValidation_t))
    names(rf_validation)<-names(map_validation_df)[-2]
    
    #make statewide valid stat from all county loocv
    
    state_statement<-as.character(paste("This",format(dataStartDate,"%B %d %Y"),"mosaic rainfall map of the State of Hawaii",
                                        "is a high spatial resolution (~250m) gridded prediction of cumulative rainfall in millimeters from midnight HST",
                                        format(dataStartDate,"%b %d %Y"),"to midnight HST",format(dataEndDate,"%b %d %Y."),
                                        #IDW statement here
                                        "This was produced by performing a climate-aided modified automatic kriging interpolations for each county extent. This kriging used log transformed daily rainfall anomaly ratios calculated from a mean monthly daily disaggrigated rf map, plus a consteint (observed mm + 1 / (mean daily mm + c) as an input. This process was done for four individually produced maps of Kauai, Honolulu (Oahu), Maui (Maui, Lanai, Molokai, & Kahoolawe) and Hawaii counties.", 
                                        "These kriging processes used", 
                                        sum(map_validation_df$stationCount), "unique station locations statewide",
                                        "and their",format(dataStartDate,"%B %d %Y"), 
                                        "recorded and/or estimated rainfall (mm) totals. Please consult each county map meta-data files for more details on map production and accuracy at the county scale.",  
                                        "A leave one out cross validation (LOOCV) of the all station data used in all four counties produced individual mean absolute error (MAE) values of:",
                                        paste(round(map_validation_df[map_validation_df$county=="KA","mae_rf_mm"],3),
                                              round(map_validation_df[map_validation_df$county=="OA","mae_rf_mm"],3),
                                              round(map_validation_df[map_validation_df$county=="MN","mae_rf_mm"],3),
                                              round(map_validation_df[map_validation_df$county=="BI","mae_rf_mm"],3),collapse=", "),
                                        "for Kauai, Honolulu (Oahu), Maui (Maui, Lanai, Molokai, & Kahoolawe) and Hawaii counties respectively.",
                                        "As a whole leave one out cross validation (LOOCV) data between observed station data and the interpolated estimate used in this map produced a mean absolute error (MAE) of:", 
                                        round(state_validation$mae_rf_mm ,3), "with a maximum observed rainfall of",map_validation_df$staRFmmMax,"mm, meaning overall this",format(dataStartDate,"%B %Y"),"statewide mosaic daily rainfall (mm) map is a",  
                                        ifelse(state_validation$qualMetric>=validRanks[1],"high quality estimate of daily rainfall.",
                                               ifelse(state_validation$qualMetric>=validRanks[2],"good quality estimate of daily rainfall.",
                                                      ifelse(state_validation$qualMetric>=validRanks[3],"moderate quality estimate of daily rainfall.",
                                                             ifelse(state_validation$qualMetric>=validRanks[4],"low quality estimate of daily rainfall, and should be used with dilligence.","lowest quality estimate of daily rainfall, and should be used with caution.")))),
                                        "All maps are subject to change as new data becomes available or unknown errors are corrected in reoccurring versions.", 
                                        "Errors in rainfall estimates do vary over space meaning any gridded rainfall value, even on higher quality maps, could still produce incorrect estimates.",
                                        "Check standard error (SE) maps to better understand spatial estimates of prediction error." 
    ))#statewide statement end
    
    keywords<-"Hawaii, Hawaiian Islands, rainfall prediction, daily precipitation, rainfall, climate, spatial interpolation, kriging"
    
    file_meta<-data.frame(dataStatement=state_statement,
                          keywords=keywords,
                          dataDate=dataStartDate,
                          dataStartDate=dataStartDate,
                          dataEndDate=dataEndDate,
                          dateProduced=Sys.Date(),
                          dataVersionType=datatype,
                          RFstationFile=filenames[1],
                          RFmmGridFile=filenames[2],
                          RFmmSEgridFile=filenames[3],
                          GeoCoordUnits="Decimal Degrees",
                          GeoCoordRefSystem=as.character(crs(grid)),
                          Xresolution=xres(grid),
                          Yresolution=yres(grid),
                          ExtentXmin=extent(grid)[1],
                          ExtentXmax=extent(grid)[2],
                          ExtentYmin=extent(grid)[3],
                          ExtentYmax=extent(grid)[4])
    
    
    final_meta<-rbind(data.frame(attribute=as.character(names(file_meta)),value=as.character(lapply(file_meta[1,], dataconvert))),
                      data.frame(attribute=paste0("Statewide_",as.character(names(state_validation)[-c(1:3)])),value=as.character(lapply(state_validation[,-c(1:3)],dataconvert))),
                      data.frame(attribute=as.character(names(rf_validation)),value=as.character(rf_validation[1,])),
                      data.frame(attribute=as.character(c("credits", "contacts")),
                                 value=as.character(c("All data produced by University of Hawaii at Manoa Water Resource Research Center (WRRC). Support provided by Change Hawaii EPSCoR funded by the National Science Foundation under EPSCoR Research Infrastructure Improvement Award #OIA-2149133",
                                                      "Matthew Lucas (mplucas@hawaii.edu), Keri Kodama (kodamak8@hawaii.edu),Ryan Longman (rlongman@hawaii.edu),Thomas Giambelluca (thomas@hawaii.edu)")))		
    )
    final_meta$attribute<-as.character(final_meta$attribute)
    return(final_meta)
  }#end statewide meta
}#end metamaker func

noRFoutputs<-function(meanRFgridwd,county,data_date,RF_day,realVals){
  message(paste(county,"No obsered rainfall."))
  
  #list to store out objects
  outList<-list()
  
  #set addC to 0
  addC<-NA #no rf
  outList[["bestC"]]<-addC
  
  #Read in Mean daily month RF grid 
  setwd(meanRFgridwd)#set mean rf raster: local pc
  rf_ras_name<-paste0(tolower(county),"_rf_mm",substr(data_date,6,7),"_daily.tif") #get mean rf raster name
  Mean_RF<- raster(rf_ras_name)#add mean rf of intended month  
  
  #Calculate anomalies and log anomalies 
  RF_day$total_rf_mm<-
    RF_day$RF_day_Anom <- 0 
  RF_day$RF_day_Anom_logK <- NA
  outList[["RF_dayBestC"]]<-RF_day
  
  #get validation metrics
  outList[["rf_validationBestC"]]<-data.frame(
    county=county,
    date=data_date,
    stationCount=nrow(RF_day),
    stationCountReal=sum(realVals),
    stationCountFill=nrow(RF_day)-sum(realVals),
    staRFmmMin=0,
    staRFmmMax=0,
    noRF=TRUE,
    addC=addC,
    fixedVario=NA,
    nugFixZero=NA,
    mod=NA,
    nugget=NA,
    range=NA,
    sill=NA,
    rsq_rf_mm=NA,
    rmse_rf_mm=NA,
    mae_rf_mm=NA,
    bias_rf_mm=NA,
    nse_rf_mm=NA,
    kge_rf_mm=NA,
    rsq_rf_anom=NA,
    rmse_rf_anom=NA,
    mae_rf_anom=NA,
    bias_rf_anom=NA,
    nse_rf_anom=NA,
    kge_rf_anom=NA
  )
  
  #make LOOCV data
  outList[["loocvBestC"]]<-data.frame(SKN=RF_day$SKN,
                                      noRF=rep(TRUE,length(RF_day$SKN)),
                                      date=rep(data_date,length(RF_day$SKN)),
                                      pred_rf=rep(0,length(RF_day$SKN)),
                                      obs_rf=RF_day$total_rf_mm,
                                      obs_anom=rep(0,length(RF_day$SKN)),
                                      pred_anom=rep(0,length(RF_day$SKN))
  )
  
  #make and save rasters
  outList[["rf_mm_ras"]]<-(Mean_RF-Mean_RF) #make 0 rf raster
  outList[["rf_mm_SE_ras"]]<-(Mean_RF/Mean_RF)*(-9999) #NA se raster
  outList[["rf_anom_ras"]]<- outList[["rf_mm_ras"]]/(Mean_RF) #rf mm anom 0
  outList[["rf_anom_SE_ras"]]<- (Mean_RF/Mean_RF)*(-9999) #NA #NA se raster
  outList[["varioModel"]]<- NA
  message(paste(county,"Map set to 0mm no krigging... complete"))
  
  #return
  return(outList)
  #end all sta rf = zero (no rainfall)
}#end no rf func

bestRFoutputs<-function(county,meanRFgridwd,data_date,zdist=0.00225,varioDFAll,RF_day,RF_day_raw,realVals){
  #list to store out objects
  outList<-list()
  
  #loops                    
  addCvec<-c(-1,0,1,10,100,1000,10000)
  dfC<-data.frame() #blank df to store all addC outputs
  allCList<-list() #blank list
  for(useVario in c(TRUE,FALSE)){ #NOTE fixed and free ON 
    for(k in 1:length(addCvec)){
      
      #define mean rf add C
      addC<-addCvec[k]
      message(paste(county,addC,"addC fix vario:",useVario,"run..."))
      
      #Read in Mean daily month RF grid 
      setwd(meanRFgridwd)#set mean rf raster: local pc
      rf_ras_name<-paste0(tolower(county),"_rf_mm",substr(data_date,6,7),"_daily.tif") #get mean rf raster name
      Mean_RF<- raster(rf_ras_name)#add mean rf of intended month  
      
      #krig no climate aid
      if(addC==-1){
        #extract Values from raster
        RF_day$RF_Mean_Extract <- raster::extract(Mean_RF, RF_day)#daily rf
        RF_day<-RF_day[!is.na(RF_day$RF_Mean_Extract),] #remove NA stations not located on mean map
        
        C<-1 #define constant C for log + C transformation
        RF_day$total_rf_mm_logC<-log(RF_day$total_rf_mm+C)
        
        if(useVario){
          #build vario
          Nugget<-varioDFAll[varioDFAll$month==month(data_date) & varioDFAll$addC==addC & varioDFAll$county==county,c("nugget")]
          Sill<-varioDFAll[varioDFAll$month==month(data_date) & varioDFAll$addC==addC & varioDFAll$county==county,c("sill")]
          Range<-varioDFAll[varioDFAll$month==month(data_date) & varioDFAll$addC==addC& varioDFAll$county==county,c("range")]
          nugFixZero<-FALSE
          vario<-try(automap::autofitVariogram(as.formula("total_rf_mm_logC ~ 1"), RF_day, model="Exp",fix.values = as.numeric(c(Nugget,Range,Sill)))) #fit mat variogram fix parameters
          if(!is.na(grep("error",vario[1],ignore.case = T)[1])){ #if error make nug 0
            Nugget<- 0 #set nugget value when needed
            vario<-automap::autofitVariogram(as.formula("total_rf_mm_logC ~ 1"), RF_day, model="Exp",fix.values = as.numeric(c(Nugget,Range,Sill))) #fit mat variogram all parameters free
            nugFixZero<-TRUE
          }
          #plot(vario)
        }else{
          nugFixZero<-FALSE
          vario<-try(automap::autofitVariogram(as.formula("total_rf_mm_logC ~ 1"), RF_day, fix.values = as.numeric(c(NA,NA,NA)))) #fit mat variogram all parameters free
          
          if(!is.na(grep("error",vario[1],ignore.case = T)[1])){ #if error make nug 0
            Nugget<- 0 #set nugget value when needed
            vario<-try(automap::autofitVariogram(as.formula("total_rf_mm_logC ~ 1"), RF_day,fix.values = as.numeric(c(Nugget,NA,NA)))) #fit mat variogram all parameters free
            nugFixZero<-TRUE
          }
          #plot(vario)
        }#end conditional free vario
        
        if(!is.na(grep("error",vario[1],ignore.case = T)[1])) next
        
        #LOOCV
        stationCount<-nrow(RF_day) #save station count with gap fill
        loocv_df<-data.frame() #blank df
        for(j in 1:nrow(RF_day)){
          check<-realFill(RF_day,RF_day_raw)[2] #check if real or gap fill val
          if(check){
            krigeLOO<-krige(as.formula("total_rf_mm_logC ~ 1") ,RF_day[-j,], RF_day[j,], model=vario$var_model)
            #back transform to rf mm and anom
            krig_logC<-krigeLOO$var1.pred
            if(!is.na(krig_logC)){
              if(krig_logC>700){krig_logC<-700} #to avoid making exp(>700) = inf
            }
            rf_mm_pred<- round((exp(krig_logC) - C),10) #back transform to get rf mm total: note not Anom
            rf_mm_obs<-RF_day[j,]$total_rf_mm
            obs_anom<-rf_mm_obs/(RF_day[j,]$RF_Mean_Extract)
            pred_anom<-rf_mm_pred/(RF_day[j,]$RF_Mean_Extract)
            #message(paste(j,"obs rf:",rf_mm_obs,"pred rf:",rf_mm_pred))
            #names(RF_day)
            cvrow<-data.frame(SKN=RF_day[j,]$SKN,
                              noRF=FALSE,
                              date=data_date,
                              pred_rf=rf_mm_pred,
                              obs_rf=rf_mm_obs,
                              obs_anom=obs_anom,
                              pred_anom=pred_anom)
            loocv_df<-rbind(loocv_df,cvrow)
          }#end real value check
        }#end loocv loop
        loocv_df[loocv_df$pred_rf<0 & !is.na(loocv_df$pred_rf),"pred_rf"]<-0 #make all ned RF pred 0
        
        #validation
        if(sum(is.na(loocv_df$pred_rf))==length(loocv_df$pred_rf)){
          rf_validation<-data.frame(county=county,date=data_date,stationCount=nrow(RF_day),stationCountReal=sum(realVals),stationCountFill=nrow(RF_day)-sum(realVals),
                                    staRFmmMin=NA,staRFmmMax=NA,noRF=NA,addC=addC,fixedVario=useVario,nugFixZero=NA,mod=NA,nugget=NA,range=NA,sill=NA,
                                    rsq_rf_mm=NA,rmse_rf_mm=NA,mae_rf_mm=NA,bias_rf_mm=NA,nse_rf_mm=NA,kge_rf_mm=NA,
                                    rsq_rf_anom=NA,rmse_rf_anom=NA,mae_rf_anom=NA,bias_rf_anom=NA,nse_rf_anom=NA,kge_rf_anom=NA
          )
          message(paste("krige produced only NA values:",county,addC,useVario))
        }else{
          rf_validation<-getValidMets(county=county,data_date=data_date,addC=addC,stationCount=stationCount,useVario=useVario,nugFixZero=nugFixZero,vario=vario,RF_day_raw=RF_day_raw,loocv_df=loocv_df)
        }#end validation
        
        dfC<-rbind(dfC,rf_validation[,c("addC","fixedVario","rsq_rf_mm","rmse_rf_mm","mae_rf_mm","bias_rf_mm","nse_rf_mm","kge_rf_mm","nugget","range","sill")])
        if(useVario){
          allCList[[k]]<-list(RF_day,vario,loocv_df,rf_validation)
          names(allCList)[k]<-paste0("add",addC,useVario)
        }else{
          allCList[[length(addCvec)+k]]<-list(RF_day,vario,loocv_df,rf_validation)
          names(allCList)[length(addCvec)+k]<-paste0("add",addC,useVario)
        }
        #end no climate aid krig
      }else{
        #climate aided krig
        Mean_RF<-(Mean_RF+addC) #ensure no zero mean rf values
        
        #extract Values from raster
        RF_day$RF_Mean_Extract <- raster::extract(Mean_RF, RF_day)#daily rf
        RF_day<-RF_day[!is.na(RF_day$RF_Mean_Extract),] #remove NA stations not located on mean map
        
        #Calculate anomalies and log anomalies 
        C<-1 #define constant C for log + C transformation
        RF_day$RF_day_Anom <- RF_day$total_rf_mm / (RF_day$RF_Mean_Extract-addC) #calc anom" rf obs/mean monthly (minus addC)
        RF_day$RF_day_Anom_logK <- log((RF_day$total_rf_mm + C * RF_day$RF_Mean_Extract)/(RF_day$RF_Mean_Extract))#Make suggested transformation a column to back transform
        
        ###Set up for interpolation###
        
        if(useVario){
          #build vario
          Nugget<-varioDFAll[varioDFAll$month==month(data_date) & varioDFAll$addC==addC & varioDFAll$county==county,c("nugget")]
          Sill<-varioDFAll[varioDFAll$month==month(data_date) & varioDFAll$addC==addC & varioDFAll$county==county,c("sill")]
          Range<-varioDFAll[varioDFAll$month==month(data_date) & varioDFAll$addC==addC& varioDFAll$county==county,c("range")]
          nugFixZero<-FALSE
          vario<-try(automap::autofitVariogram(as.formula("RF_day_Anom_logK ~ 1"), RF_day, model="Exp",fix.values = as.numeric(c(Nugget,Range,Sill)))) #fit mat variogram fix parameters
          if(!is.na(grep("error",vario[1],ignore.case = T)[1])){ #if error make nug 0
            Nugget<- 0 #set nugget value when needed
            vario<-automap::autofitVariogram(as.formula("RF_day_Anom_logK ~ 1"), RF_day, model="Exp",fix.values = as.numeric(c(Nugget,Range,Sill))) #fit mat variogram all parameters free
            nugFixZero<-TRUE
          }
          #plot(vario)
        }else{
          nugFixZero<-FALSE
          vario<-try(automap::autofitVariogram(as.formula("RF_day_Anom_logK ~ 1"), RF_day, fix.values = as.numeric(c(NA,NA,NA)))) #fit mat variogram all parameters free
          if(!is.na(grep("error",vario[1],ignore.case = T)[1])){ #if error make nug 0
            Nugget<- 0 #set nugget value when needed           
            nugFixZero<-TRUE
            vario<-try(automap::autofitVariogram(as.formula("RF_day_Anom_logK ~ 1"), RF_day, fix.values = as.numeric(c(Nugget,NA,NA))))#fit mat variogram all parameters free
          }
          #plot(vario)
        }#end conditional free vario
        if(!is.na(grep("error",vario[1],ignore.case = T)[1])) next  #if still error skip rest of loop iteration for county
        
        ## make loocv df
        stationCount<-nrow(RF_day) #save station count with gap fill
        loocv_df<-data.frame() #blank df
        for(j in 1:nrow(RF_day)){
          check<-realFill(RF_day,RF_day_raw)[j] #check if real or gap fill val
          if(check){
            krigeLOO<-krige(as.formula("RF_day_Anom_logK ~ 1") ,RF_day[-j,], RF_day[j,], model=vario$var_model)
            #back transform to rf mm and anom
            krig_logC_anom<-krigeLOO$var1.pred
            if(!is.na(krig_logC_anom)){
              if(krig_logC_anom>700){krig_logC_anom<-700} #to avoid making exp(>700) = inf
            }
            rf_mm_pred<- round((exp(krig_logC_anom) - C)*RF_day[j,]$RF_Mean_Extract,10) #back transform to get rf mm total: note not Anom
            rf_mm_obs<-RF_day[j,]$total_rf_mm
            obs_anom<-RF_day[j,]$RF_day_Anom
            pred_anom<-rf_mm_pred/(RF_day[j,]$RF_Mean_Extract-addC)
            #message(paste(j,"obs rf:",rf_mm_obs,"pred rf:",rf_mm_pred))
            #names(RF_day)
            cvrow<-data.frame(SKN=RF_day[j,]$SKN,
                              noRF=FALSE,
                              date=data_date,
                              pred_rf=rf_mm_pred,
                              obs_rf=rf_mm_obs,
                              obs_anom=obs_anom,
                              pred_anom=pred_anom)
            loocv_df<-rbind(loocv_df,cvrow)
          }#end real value check
        }#end loocv loop
        
        loocv_df[loocv_df$pred_rf<0 & !is.na(loocv_df$pred_rf),"pred_rf"]<-0 #make all neg RF pred 0
        
        #validation metrics
        if(sum(is.na(loocv_df$pred_rf))==length(loocv_df$pred_rf)){
          rf_validation<-data.frame(county=county,date=data_date,stationCount=nrow(RF_day),stationCountReal=sum(realVals),stationCountFill=nrow(RF_day)-sum(realVals),
                                    staRFmmMin=NA,staRFmmMax=NA,noRF=NA,addC=addC,fixedVario=useVario,nugFixZero=NA,mod=NA,nugget=NA,range=NA,sill=NA,
                                    rsq_rf_mm=NA,rmse_rf_mm=NA,mae_rf_mm=NA,bias_rf_mm=NA,nse_rf_mm=NA,kge_rf_mm=NA,
                                    rsq_rf_anom=NA,rmse_rf_anom=NA,mae_rf_anom=NA,bias_rf_anom=NA,nse_rf_anom=NA,kge_rf_anom=NA
          )
          message(paste("krige produced only NA values:",county,addC,useVario))
        }else{
          rf_validation<-getValidMets(county=county,data_date=data_date,addC=addC,stationCount=stationCount,useVario=useVario,nugFixZero=nugFixZero,vario=vario,RF_day_raw=RF_day_raw,loocv_df=loocv_df)
        }#end validation
        
        #gather addC items
        dfC<-rbind(dfC,rf_validation[,c("addC","fixedVario","rsq_rf_mm","rmse_rf_mm","mae_rf_mm","bias_rf_mm","nse_rf_mm","kge_rf_mm","nugget","range","sill")])
        if(useVario){
          allCList[[k]]<-list(RF_day,vario,loocv_df,rf_validation)
          names(allCList)[k]<-paste0("add",addC,useVario)
        }else{
          allCList[[length(addCvec)+k]]<-list(RF_day,vario,loocv_df,rf_validation)
          names(allCList)[length(addCvec)+k]<-paste0("add",addC,useVario)
        }#gather c runs end
      }#end climate aided krig
    }#end addC loop
  }#useVario T/F loop
  
  ##choose best C run ranked best krig
  dfC$listName<-paste0("add",dfC$addC,dfC$fixedVario)
  dfC$zdist<-zdist
  
  #if free vario subset runs with nug/sill ratio and small ranges
  nugsill<-(dfC$nugget/dfC$sill)<0.3
  rangeLim<-dfC$range>0.04
  
  dfCFree<-dfC[nugsill & rangeLim & !dfC$fixedVario,] #subset free
  dfCFix<-dfC[dfC$fixedVario,] #subset fixed
  dfC<-rbind(dfCFix,dfCFree) #redefine dfc
  
  #conditional has data 
  if(sum(!is.na(dfC$mae_rf_mm))>0 & nrow(dfC)>0){
    
    # #get bestC meanScore
    # bestRun<-dfC[which.max(dfC$meanScore),"listName"]
    # bestC<-dfC[which.max(dfC$meanScore),"addC"]
    # bestFixVario<-dfC[which.max(dfC$meanScore),"fixedVario"]
    # bestCList<-allCList[[bestRun]]
    
    #get bestC rmse
    bestRun<-dfC[which.min(dfC$mae_rf_mm),"listName"]
    bestC<-dfC[which.min(dfC$mae_rf_mm),"addC"]
    bestFixVario<-dfC[which.min(dfC$mae_rf_mm),"fixedVario"]
    bestCList<-allCList[[bestRun]]
    
    #get bestC objects from best cc list 
    RF_dayBestC<-bestCList[[1]]
    varioBestC<-bestCList[[2]]
    loocvBestC<-bestCList[[3]]
    rf_validationBestC<-bestCList[[4]]
    message(county," best C run: ",bestRun," selected")
    
    ###krige best C grid data ###
    #re-read in Mean daily month RF grid 
    setwd(meanRFgridwd)#set mean rf raster: local pc
    rf_ras_name<-paste0(tolower(county),"_rf_mm",substr(data_date,6,7),"_daily.tif") #get mean rf raster name
    Mean_RF<- raster(rf_ras_name)#add mean rf of intended month  
    
    #kriging
    temppoints<-SpatialPoints(as.data.frame(Mean_RF,xy=T,na.rm=T)[,c(1,2)]) #make spatial points data frame with x y coords from mask and remove NA pixels
    if(bestC==-1){ #make non-climate aid raster
      
      #raw krige
      krigeObj<-krige(as.formula("total_rf_mm_logC ~ 1") ,RF_dayBestC, temppoints, model=varioBestC$var_model)
      krig_logC<-rasterize(krigeObj, Mean_RF, krigeObj$var1.pred) #make krig log points into raster
      krig_logC_SE<-rasterize(krigeObj, Mean_RF, krigeObj$var1.var) #make krig SD logk points into raster
      
      #convert to rf mm units
      rf_mm_ras<- round((exp(krig_logC) - C),10) #back transform to get rf mm total: note not Anom
      rf_mm_SE_ras<- round((exp(krig_logC_SE) - C)*Mean_RF,10) #back transform SE to get rf mm total: note not Anom
      rf_mm_ras[rf_mm_ras<0]<-0 #edit neg values to zero
      
      #make anom rasters
      rf_anom_ras<- rf_mm_ras/Mean_RF #back transform to get rf mm anom
      rf_anom_SE_ras<- rf_mm_SE_ras/Mean_RF #back transform to get anom SE
      
      #raw krig end
    }else{ #make clim aid raster
      
      Mean_RF<-(Mean_RF+bestC) #add best addC to mean
      
      #clim aid kriging
      krigeObj<-krige(as.formula("RF_day_Anom_logK ~ 1") ,RF_dayBestC, temppoints, model=varioBestC$var_model)
      krig_logC_anom<-rasterize(krigeObj, Mean_RF, krigeObj$var1.pred) #make krig log points into raster
      krig_logC_anom_SE<-rasterize(krigeObj, Mean_RF, krigeObj$var1.var) #make krig SD logk points into raster
      krig_logC_anom[krig_logC_anom>=708]<-700
      #plot(krig_logC_anom)
      
      #convert to rf mm units
      rf_mm_ras<- round((exp(krig_logC_anom) - C)*Mean_RF,10) #back transform to get rf mm total: note not Anom
      rf_mm_SE_ras<- round((exp(krig_logC_anom_SE) - C)*Mean_RF,10) #back transform SE to get rf mm total: note not Anom
      rf_mm_ras[rf_mm_ras<0]<-0 #edit neg values to zero
      
      #make anom rasters
      rf_anom_ras<- rf_mm_ras/(Mean_RF-bestC) #back transform to get rf mm anom
      rf_anom_SE_ras<- rf_mm_SE_ras/(Mean_RF-bestC) #back transform to get anom SE
    }#clim aid raster end
    message("rasters made ",county)
    
    #function output
    outList[["bestC"]]<-bestC
    outList[["RF_dayBestC"]]<-RF_dayBestC
    outList[["rf_validationBestC"]]<-rf_validationBestC
    outList[["loocvBestC"]]<-loocvBestC
    outList[["rf_mm_ras"]]<-rf_mm_ras
    outList[["rf_mm_SE_ras"]]<-rf_mm_SE_ras
    outList[["rf_anom_ras"]]<-rf_anom_ras
    outList[["rf_anom_SE_ras"]]<-rf_anom_SE_ras
    outList[["varioBestC"]]<-varioBestC
  }else{ 
    #no krige 
    outList[["bestC"]]<-NA
    #print(dfC)
    stop(paste("null krige",data_date,county)) #stop function kick error
    #try IDW here
  }    
  return(outList)
}#end obs rf best krig

dailyRFkrig<-function(rfdailyDFmaster,rfdailyRawDFmaster,varioDFAll,data_date,meanRFgridwd,outdir,dataVersion,testRun=NA){
  
  require(automap)
  require(sp)
  require(raster)
  require(gstat)
  require(hydroGOF)
  require(Metrics)
  require(lubridate)
  
  rgdal::setCPLConfigOption("GDAL_PAM_ENABLED", "FALSE") #set rgdal config to not writ aux.xml file
  
  #define loop counties
  counties<-c("BI","MN","OA","KA")
  
  #define data date as date
  data_date<-base::as.Date(data_date)
  
  #get date cols from rf df
  date_cols<-grep("X",names(rfdailyDFmaster))
  date_cols_raw<-grep("X",names(rfdailyRawDFmaster))
  
  #get meta cols
  metacols<-1:(min(date_cols)-1)
  
  t1<-Sys.time()
  domRun<-0 #keep track of run
  for(dom in 1:length(counties)){
    domRun<-domRun+1 #relative dom run
    #load rf station data again
    rfdailyDF<-rfdailyDFmaster
    rfdailyDFRaw<-rfdailyRawDFmaster
    
    #define county from islands and subset
    if(!is.na(testRun)){
      county=testRun #test run county define
    }else{
      county<-counties[dom] #define co for loop
    }
    if(county=="MN"){islands<-c("MA","LA","MO","KO")}else{islands<-county}
    rfdailyDF<-rfdailyDF[rfdailyDF$Island %in% islands,] #subset county data
    rfdailyDFRaw<-rfdailyDFRaw[rfdailyDFRaw$Island %in% islands,] #subset county data
    
    #subset day and rename rf date col
    dateCol<-date_cols[which(format(data_date,"X%Y.%m.%d")==names(rfdailyDF)[date_cols])] #define date column from data_date and date_cols
    dateColRaw<-date_cols_raw[which(format(data_date,"X%Y.%m.%d")==names(rfdailyDFRaw)[date_cols_raw])] #define date column from data_date and date_cols
    RF_day<-rfdailyDF[!is.na(rfdailyDF[,dateCol]),c(metacols,dateCol)] #subset meta and dateCol stations with rf data on the day (not na)
    RF_day_raw<-rfdailyDFRaw[!is.na(rfdailyDFRaw[,dateColRaw]),c(1,dateColRaw)] #subset SKNs with rf data on data date
    names(RF_day)[max(metacols)+1]<-"total_rf_mm" #rename date col as total_rf_mm
    names(RF_day_raw)[2]<-"total_rf_mm" #rename date col as total_rf_mm
    realVals<-realFill(RF_day,RF_day_raw) #real daily rf values only subset vec
    #print(RF_day) #check final rf day df
    
    #make spatial data frame with Lat and Lon
    RF_day$x<-RF_day$LON
    RF_day$y<-RF_day$LAT
    coordinates(RF_day) <- ~x+y
    crs(RF_day)<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" #add crs 
    
    #check station count is greater then N
    if(nrow(RF_day)<3){ #is there less then N stations?
      warning(paste(nrow(RF_day),"is not enough stations to krige",data_date,county))
      rf_validation<-data.frame(
        county=county,
        date=data_date,
        stationCount=nrow(RF_day),
        stationCountReal=sum(realVals),
        stationCountFill=nrow(RF_day)-sum(realVals),
        staRFmmMin=NA,
        staRFmmMax=NA,
        noRF=NA,
        addC=addC,
        fixedVario=NA,
        nugFixZero=NA,
        mod=NA,
        nugget=NA,
        range=NA,
        sill=NA,
        rsq_rf_mm=NA,
        rmse_rf_mm=NA,
        mae_rf_mm=NA,
        bias_rf_mm=NA,
        nse_rf_mm=NA,
        kge_rf_mm=NA,
        rsq_rf_anom=NA,
        rmse_rf_anom=NA,
        mae_rf_anom=NA,
        bias_rf_anom=NA,
        nse_rf_anom=NA,
        kge_rf_anom=NA
      )
      return(rf_validation)
    }else{ #more then 3 station conditional
      
      #set up addC loop, fix vario condition and while loop condition  
      zdist <- 0.00225 #set zero dist
      m=1
      
      #remove dup location
      if(m==1){
        net_priority<-c("Mesonet","HaleNet","CraterNet","Little,HaleNet","HavoNet","HIPPNET","USCRN","ESRL/GMD","NWS","USGS","HydroNet-UaNet","SCAN","NREL","RAWS","STATE","NREM","PrivateObs","HC&S","COOP","CoCoRaHS","CoCoRaHs")
        RF_day_real<-RF_day[realVals,]#sub real vals
        RF_day_fill<-RF_day[!realVals,]#sub fill vals
        RF_day_real<-RF_day_real[order(match(RF_day_real$Network, net_priority)),]  #reorder by network for real values
        RF_day_fill<-RF_day_fill[order(match(RF_day_fill$Network, net_priority)),]  #reorder by network for fill values
        RF_day<-rbind(RF_day_real,RF_day_fill) #combine together real and fill
      }
      crs(RF_day)<-NA #remove CRS for vario and location
      dupLocal<-zerodist(RF_day,zero = zdist*m) #distance within zdist * m
      if(nrow(dupLocal)>0){RF_day <- RF_day[-dupLocal[,2],]} #remove dup local
      
      #recheck station min/max
      staRFmin<-min(RF_day$total_rf_mm) #save station min stat
      staRFmax<-max(RF_day$total_rf_mm) #save station max stat
      staRFmaxNoFill<-max(RF_day_real$total_rf_mm) #save station max no fill station stat
      
      if(staRFmaxNoFill==0 | staRFmax==0){
        #no RF function and outputs
        message("no rf observed...")
        bestCList<-noRFoutputs(meanRFgridwd,county,data_date,RF_day,realVals)
      }else{ 
        #station observed rainfall >0 else
        message("rf observed attempting best krig...")
        bestCList<-bestRFoutputs(county,meanRFgridwd,data_date,zdist,varioDFAll,RF_day,RF_day_raw,realVals)
      }#end RF=0 if else
      
      #get bestC objects from best cc list 
      bestC<-bestCList[["bestC"]]
      RF_dayBestC<-bestCList[["RF_dayBestC"]]
      varioBestC<-bestCList[["varioBestC"]]
      loocvBestC<-bestCList[["loocvBestC"]]
      rf_validationBestC<-bestCList[["rf_validationBestC"]]
      rf_mm_ras<-bestCList[["rf_mm_ras"]]
      rf_mm_SE_ras<-bestCList[["rf_mm_SE_ras"]]
      rf_anom_ras<-bestCList[["rf_anom_ras"]]
      rf_anom_SE_ras<-bestCList[["rf_anom_SE_ras"]]
      bestFixVario<-rf_validationBestC$fixedVario
      
      ##county plots##
      dir.create(paste0(outdir,"/plots"), showWarnings = FALSE)
      dir.create(paste0(outdir,"/plots/daily"), showWarnings = FALSE)
      
      #variogram
      dir.create(paste0(outdir,"/plots/daily/variogram"), showWarnings = FALSE)
      dir.create(paste0(outdir,"/plots/daily/variogram/county"), showWarnings = FALSE)
      dir.create(paste0(outdir,"/plots/daily/variogram/county/",county), showWarnings = FALSE)
      setwd(paste0(outdir,"/plots/daily/variogram/county/",county))
      if(!rf_validationBestC$noRF){
        subVG<-paste0(county," ",format(data_date,"%Y-%m-%d"),"; Fix Vario:",bestFixVario,"; BestC:",bestC,"; RSQ:",round(rf_validationBestC$rsq_rf_mm,2),"; RMSE:",round(rf_validationBestC$rmse_rf_mm,2),"; MAE:",round(rf_validationBestC$mae_rf_mm,2),"; BIAS:",round(rf_validationBestC$bias_rf_mm,2))
        bitmap(file = paste0(county,"_vario_",format(data_date,"%Y%m%d"),".jpg"),width=7,height=5,units="in",res=300,type="jpeg")
          plot(varioBestC,sub=subVG)
          dev.off()
        message(paste(county,"vario plotted!"))
      }
      
      #rf mm plotting
      rf_mm_ras[rf_mm_ras<0]<-0 #edit neg values to zero
      #rf_mm_ras_zero<-rf_mm_ras<=0.001
      #rf_mm_ras_zero[rf_mm_ras_zero == 0] <- NA
      mapMax<-max(c(max(RF_dayBestC$total_rf_mm),maxValue(rf_mm_ras)))
      pntCol <- rainbow(100,end=0.8)[as.numeric(cut(c(0,mapMax,RF_dayBestC$total_rf_mm),breaks = 100))[-c(1,2)]]
      subText<-paste0("Fix Vario:",bestFixVario,"; BestC:",bestC,";  RSQ:",round(rf_validationBestC$rsq_rf_mm,2),";  RMSE:",round(rf_validationBestC$rmse_rf_mm,2),";  MAE:",round(rf_validationBestC$mae_rf_mm,2),";  BIAS:",round(rf_validationBestC$bias_rf_mm,2))
      #rf_mm_poly_zero<-rasterToPolygons(rf_mm_ras_zero, digits=8, dissolve=TRUE)
      
      #save rf plot
      # dir.create(paste0(outdir,"/plots"), showWarnings = FALSE)
      # dir.create(paste0(outdir,"/plots/daily"), showWarnings = FALSE)
      dir.create(paste0(outdir,"/plots/daily/rf_mm"), showWarnings = FALSE)
      dir.create(paste0(outdir,"/plots/daily/rf_mm/county"), showWarnings = FALSE)
      dir.create(paste0(outdir,"/plots/daily/rf_mm/county/",county), showWarnings = FALSE)
      setwd(paste0(outdir,"/plots/daily/rf_mm/county/",county))
      bitmap(file = paste0(county,"_rf_mm_",format(data_date,"%Y%m%d"),".jpg"),width=7,height=5,units="in",res=300,type="jpeg")
      plot(rf_mm_ras,col=rainbow(100,end=0.8),zlim=c(0,mapMax),main=paste("Daily RF mm:",data_date),sub=subText)
      points(RF_dayBestC[RF_dayBestC$SKN %in% RF_day_raw$SKN,], cex=1.25, pch=21, bg=pntCol[RF_dayBestC$SKN %in% RF_day_raw$SKN],col="black")
      points(RF_dayBestC[!RF_dayBestC$SKN %in% RF_day_raw$SKN,], cex=1.25, pch=22, bg=pntCol[!RF_dayBestC$SKN %in% RF_day_raw$SKN],col="black")
      points(RF_dayBestC[RF_dayBestC$total_rf_mm==0,], cex=.75, pch=4, col="black")
      dev.off()
      
      #save anom plot
      subTextAnom<-paste0("Fix Vario:",bestFixVario,"; BestC:",bestC,";  RSQ:",round(rf_validationBestC$rsq_rf_anom,2),";  RMSE:",round(rf_validationBestC$rmse_rf_anom,2),";  MAE:",round(rf_validationBestC$mae_rf_anom,2),";  BIAS:",round(rf_validationBestC$bias_rf_anom,2))
      dir.create(paste0(outdir,"/plots/daily/anom"), showWarnings = FALSE)
      dir.create(paste0(outdir,"/plots/daily/anom/county"), showWarnings = FALSE)
      dir.create(paste0(outdir,"/plots/daily/anom/county/",county), showWarnings = FALSE)
      setwd(paste0(outdir,"/plots/daily/anom/county/",county))
      bitmap(file = paste0(county,"_anom_",format(data_date,"%Y%m%d"),".jpg"),width=7,height=5,units="in",res=300,type="jpeg")
      plot(rf_anom_ras,col=rev(topo.colors(100)),main=paste("Daily RF Anomaly:",data_date),sub=subTextAnom)
      dev.off()
      
      message("plots made and saved ",county)
      #end plots
      
      ##write and store rf mm raster
      dir.create(paste0(outdir,"/tiffs"), showWarnings = FALSE)
      dir.create(paste0(outdir,"/tiffs/daily"), showWarnings = FALSE)
      dir.create(paste0(outdir,"/tiffs/daily/county"), showWarnings = FALSE)
      dir.create(paste0(outdir,"/tiffs/daily/county/rf_mm"), showWarnings = FALSE)
      dir.create(paste0(outdir,"/tiffs/daily/county/rf_mm/",county), showWarnings = FALSE)
      setwd(paste0(outdir,"/tiffs/daily/county/rf_mm/",county))
      rfmmfilename<-paste0(format(data_date,"%Y%m%d_"),tolower(county),"_rf_mm.tif")
      writeRaster(rf_mm_ras,rfmmfilename,overwrite=TRUE)
      if(domRun==1){state_rf_mm_ras<-rf_mm_ras}else{state_rf_mm_ras<-mosaic(state_rf_mm_ras,rf_mm_ras,fun=max)} #build state mosiac raster
      
      #write and store rf mm SE raster
      dir.create(paste0(outdir,"/tiffs/daily/county/rf_mm_se"), showWarnings = FALSE)
      dir.create(paste0(outdir,"/tiffs/daily/county/rf_mm_se/",county), showWarnings = FALSE)
      setwd(paste0(outdir,"/tiffs/daily/county/rf_mm_se/",county))
      rfmmSEfilename<-paste0(format(data_date,"%Y%m%d_"),tolower(county),"_rf_mm_SE.tif")
      writeRaster(rf_mm_SE_ras,rfmmSEfilename,overwrite=TRUE)
      if(domRun==1){state_rf_mm_SE_ras<-rf_mm_SE_ras}else{state_rf_mm_SE_ras<-mosaic(state_rf_mm_SE_ras,rf_mm_SE_ras,fun=max)} #build state mosiac raster
      
      #write and store anom raster
      dir.create(paste0(outdir,"/tiffs/daily/county/anom"), showWarnings = FALSE)
      dir.create(paste0(outdir,"/tiffs/daily/county/anom/",county), showWarnings = FALSE)
      setwd(paste0(outdir,"/tiffs/daily/county/anom/",county))
      rfAnomfilename<-paste0(format(data_date,"%Y%m%d_"),tolower(county),"_anom.tif")
      writeRaster(rf_anom_ras,rfAnomfilename,overwrite=TRUE)
      if(domRun==1){state_rf_anom_ras<-rf_anom_ras}else{state_rf_anom_ras<-mosaic(state_rf_anom_ras,rf_anom_ras,fun=max)} #build state mosiac raster
      
      #write and store anom SE raster
      dir.create(paste0(outdir,"/tiffs/daily/county/anom_se"), showWarnings = FALSE)
      dir.create(paste0(outdir,"/tiffs/daily/county/anom_se/",county), showWarnings = FALSE)
      setwd(paste0(outdir,"/tiffs/daily/county/anom_se/",county))
      rfAnomSEfilename<-paste0(format(data_date,"%Y%m%d_"),tolower(county),"_anom_SE.tif")
      writeRaster(rf_anom_SE_ras,rfAnomSEfilename,overwrite=TRUE)
      if(domRun==1){state_rf_anom_SE_ras<-rf_anom_SE_ras}else{state_rf_anom_SE_ras<-mosaic(state_rf_anom_SE_ras,rf_anom_SE_ras,fun=max)} #build state mosiac raster
      
      message("rasters saved ",county)
      
      ##save anom krig input
      dir.create(paste0(outdir,"/tables"), showWarnings = FALSE)
      dir.create(paste0(outdir,"/tables/station_data"), showWarnings = FALSE)
      dir.create(paste0(outdir,"/tables/station_data/daily"), showWarnings = FALSE)
      dir.create(paste0(outdir,"/tables/station_data/daily/krigInput"), showWarnings = FALSE)
      dir.create(paste0(outdir,"/tables/station_data/daily/krigInput/county"), showWarnings = FALSE)
      dir.create(paste0(outdir,"/tables/station_data/daily/krigInput/county/",county), showWarnings = FALSE)
      setwd(paste0(outdir,"/tables/station_data/daily/krigInput/county/",county))
      krigInFilename<-paste0(format(data_date,"%Y%m%d_"),tolower(county),"_rf_krig_input.csv")
      write.csv(as.data.frame(RF_dayBestC),krigInFilename,row.names=F)
      
      if(domRun==1){stateKrigInput<-as.data.frame(RF_dayBestC)}else{stateKrigInput<-rbindAll(stateKrigInput,as.data.frame(RF_dayBestC))}
      message("krig input saved ",county)
      
      ##save loocv
      dir.create(paste0(outdir,"/tables/validation"),showWarnings = FALSE)
      dir.create(paste0(outdir,"/tables/validation/daily"),showWarnings = FALSE)
      dir.create(paste0(outdir,"/tables/validation/daily/loocv"),showWarnings = FALSE)
      dir.create(paste0(outdir,"/tables/validation/daily/loocv/county"), showWarnings = FALSE)
      dir.create(paste0(outdir,"/tables/validation/daily/loocv/county/",county), showWarnings = FALSE)
      setwd(paste0(outdir,"/tables/validation/daily/loocv/county/",county))
      loocvFilename<-paste0(format(data_date,"%Y%m%d_"),tolower(county),"_rf_loocv.csv")
      write.csv(loocvBestC,loocvFilename,row.names=F)
      if(domRun==1){stateloocv<-loocvBestC}else{stateloocv<-rbindAll(stateloocv,loocvBestC)}
      
      message("loocv saved ",county)
      
      ##make and save metadata
      filenames<-c(krigInFilename,rfmmfilename,rfmmSEfilename)
      countyMeta<-metamaker(map_validation_df=rf_validationBestC,
                            grid=rf_mm_ras,
                            filenames=filenames,
                            datatype=dataVersion,
                            rfDay=data_date,
                            statewide=F)
      
      #save metadata
      dir.create(paste0(outdir,"/tables/metadata"), showWarnings = FALSE)
      dir.create(paste0(outdir,"/tables/metadata/daily"), showWarnings = FALSE)
      dir.create(paste0(outdir,"/tables/metadata/daily/county"), showWarnings = FALSE)
      dir.create(paste0(outdir,"/tables/metadata/daily/county/",county), showWarnings = FALSE)
      setwd(paste0(outdir,"/tables/metadata/daily/county/",county))
      coMetaFilename<-paste0(format(data_date,"%Y%m%d_"),tolower(county),"_rf_meta.txt")
      write.table(countyMeta, coMetaFilename,sep ="\t", row.names = F, quote = F)
      
      message("metadata made and saved ",county)
      
      #save validation row
      if(domRun==1){allCounties_rf_validation<-rf_validationBestC}else{allCounties_rf_validation<-rbind(allCounties_rf_validation,rf_validationBestC)}
      
      message(county," complete")
    }#less then 3 station else end
    if(!is.na(testRun)){
      break # break loop if test ie: only run test county
      stop("county test run ended")
    }#end test county run
  }#end county loop
  t2<-Sys.time()
  t2-t1
  
  #save STATEWIDE files
  message("saving statewide files...")
  
  #all counties validation
  dir.create(paste0(outdir,"/tables/validation"), showWarnings = FALSE)
  dir.create(paste0(outdir,"/tables/validation/daily"), showWarnings = FALSE)
  dir.create(paste0(outdir,"/tables/validation/daily/validate"), showWarnings = FALSE)
  dir.create(paste0(outdir,"/tables/validation/daily/validate"), showWarnings = FALSE)
  dir.create(paste0(outdir,"/tables/validation/daily/validate/counties"), showWarnings = FALSE)
  setwd(paste0(outdir,"/tables/validation/daily/validate/counties"))
  allCountiesvalidFilename<-paste0(format(data_date,"%Y%m%d"),"_Counties_rf_daily_validation.csv")
  write.csv(allCounties_rf_validation,allCountiesvalidFilename,row.names=F)
  
  # #reorder date
  # allCounties_rf_validation<-read.csv(allCountiesvalidFilename)
  # allCounties_rf_validation$date<-as.Date(allCounties_rf_validation$date)
  # allCounties_rf_validation <- allCounties_rf_validation[order(allCounties_rf_validation$date),]
  # write.csv(allCounties_rf_validation,allCountiesvalidFilename,row.names=F)
  
  #state krige input
  dir.create(paste0(outdir,"/tables/station_data/daily/krigInput/statewide"), showWarnings = FALSE)
  setwd(paste0(outdir,"/tables/station_data/daily/krigInput/statewide"))
  statekrigInFilename<-paste0(format(data_date,"%Y%m%d_"),"Statewide_rf_krig_input.csv")
  write.csv(stateKrigInput,statekrigInFilename,row.names=F)
  
  #state loocv
  dir.create(paste0(outdir,"/tables/validation/daily/loocv"), showWarnings = FALSE)
  dir.create(paste0(outdir,"/tables/validation/daily/loocv/statewide"), showWarnings = FALSE)
  setwd(paste0(outdir,"/tables/validation/daily/loocv/statewide"))
  stateloocvFilename<-paste0(format(data_date,"%Y%m%d"),"_Statewide_rf_loocv.csv")
  write.csv(stateloocv,stateloocvFilename,row.names=F)
  
  #make statewide validation table row to append
  statewide_rf_validation<-getValidMets(statewide = T,data_date=data_date,stationCount=nrow(stateKrigInput),loocv_df=stateloocv)
  dir.create(paste0(outdir,"/tables/validation/daily/validate/statewide"), showWarnings = FALSE)
  setwd(paste0(outdir,"/tables/validation/daily/validate/statewide"))
  statevalidFilename<-paste0(format(data_date,"%Y%m%d_"),"Statewide_rf_daily_validation.csv")
  write.csv(statewide_rf_validation,statevalidFilename,row.names=F)
  
  message("statewide tables saved")
  
  #STATEWIDE RASTERS
  #write statewide RF raster files (not anoms)
  dir.create(paste0(outdir,"/tiffs/daily/statewide"), showWarnings = FALSE)
  dir.create(paste0(outdir,"/tiffs/daily/statewide/rf_mm"), showWarnings = FALSE)
  setwd(paste0(outdir,"/tiffs/daily/statewide/rf_mm"))
  staterfmmfilename<-paste0(format(data_date,"%Y%m%d_"),"Statewide_rf_mm.tif")
  writeRaster(state_rf_mm_ras,staterfmmfilename,overwrite=TRUE)
  
  #rf SE
  dir.create(paste0(outdir,"/tiffs/daily/statewide/rf_mm_se"), showWarnings = FALSE)
  setwd(paste0(outdir,"/tiffs/daily/statewide/rf_mm_se"))
  staterfmmSEfilename<-paste0(format(data_date,"%Y%m%d_"),"Statewide_rf_mm_SE.tif")
  writeRaster(state_rf_mm_SE_ras,staterfmmSEfilename,overwrite=TRUE)
  
  message("statewide rf rasters saved")
  
  
  ## Statewide rf map
  subText<-paste0("RSQ:",round(statewide_rf_validation$rsq_rf_mm,2),";   RMSE:",round(statewide_rf_validation$rmse_rf_mm,2),";   MAE:",round(statewide_rf_validation$mae_rf_mm,2),";   BIAS:",round(statewide_rf_validation$bias_rf_mm,2))
  
  #save state rf plot
  dir.create(paste0(outdir,"/plots/daily/rf_mm/statewide"), showWarnings = FALSE)
  setwd(paste0(outdir,"/plots/daily/rf_mm/statewide"))
  bitmap(file = paste0("Statewide_rf_mm_",format(data_date,"%Y%m%d"),".jpg"),width=7,height=5,units="in",res=300,type="jpeg")
  plot(state_rf_mm_ras,col=rainbow(100,end=0.8),main=paste("Daily RF mm:",data_date),sub=subText)
  dev.off()
  message("Statewide rf map saved")
  
  # #state rf set legend plot
  # dir.create(paste0(outdir,"/plots/daily/rf_mm/statewide_set_leg"), showWarnings = FALSE)
  # state_rf_mm_ras2<-state_rf_mm_ras
  # state_rf_mm_ras2[state_rf_mm_ras2 >= 75] <- 75 #trunc upper mm for plot
  # setwd(paste0(outdir,"/plots/daily/rf_mm/statewide_AGU"))
  # bitmap(file = paste0("Statewide_rf_100mm_",format(data_date,"%Y%m%d"),".jpg"),width=1920,height=1080,units="px",res=300,type="jpeg")
  # plot(state_rf_mm_ras2,col=rainbow(100,end=0.8),main=paste("Daily RF mm:",data_date),sub=subText,zlim=c(0,75),legend.args = list(text = ' >75', side = 3))
  # plot(state_rf_mm_ras2,col=rainbow(100,end=0.8), zlim=c(0,75),legend.only=TRUE, legend.args=list(text='Rainfall (mm)', side=4, font=2, line=2.5, cex=0.8))
  # dev.off()
  # message("Statewide rf agu map saved")
  
  #state metadata
  statefilenames<-c(statekrigInFilename,staterfmmfilename,staterfmmSEfilename,stateloocvFilename)
  
  StateMeta<-metamaker(map_validation_df=allCounties_rf_validation,
                       state_validation=statewide_rf_validation,
                       grid=rf_mm_ras,
                       filenames=statefilenames,
                       datatype="prelimanary",
                       rfDay=data_date,
                       statewide=TRUE,
                       loocv_df =stateloocv)
  #save metadata
  dir.create(paste0(outdir,"/tables/metadata/daily/statewide"), showWarnings = FALSE)
  setwd(paste0(outdir,"/tables/metadata/daily/statewide"))
  StateMetaFilename<-paste0(format(data_date,"%Y%m%d_"),"Statewide_rf_meta.txt")
  write.table(StateMeta, StateMetaFilename, sep ="\t", row.names = F, quote = F)
  message("Statewide metadata made and saved")
  
  #function return
  return(statewide_rf_validation)
}


