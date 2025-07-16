rm(list = ls())



#install.packages(c("colorRamps","dplyr","randomForest","foreach","doParallel","svMisc","Metrics","gstat","automap","raster","rgdal","rgeos","lubridate"))

require(rgeos)

require(rgdal)

require(raster)

require(automap)

require(gstat)

require(Metrics)

require(svMisc)

require(doParallel)

require(foreach)

require(randomForest)

require(dplyr)

require(colorRamps)

require(lubridate)



#pkg vec for para loop

pkgs_vec<-c("colorRamps","dplyr","randomForest","foreach","doParallel","svMisc","Metrics","gstat","automap","raster","rgdal","rgeos","lubridate")



#dirs

masterDir<-"/home/hawaii_climate_products_container/preliminary/rainfall" #master dir 

runVersionDir<-paste0(masterDir,"/data_outputs/") 

dependDir<-paste0(masterDir,"/dependencies/monthly")

codeDir<-paste0(masterDir,"/code/monthly")

inDir<-paste0(masterDir,"/data_outputs/tables/station_data/monthly/partial_filled/statewide")

krigInputOut<-paste0(masterDir,"/data_outputs/tables/station_data/monthly/krigInput/statewide")

krigInputOutCo<-paste0(masterDir,"/data_outputs/tables/station_data/monthly/krigInput/county")

plotsOut<-paste0(masterDir,"/data_outputs/plots/monthly")

plotsMultiOut<-paste0(masterDir,"/data_outputs/plots/monthly/multiplot")

plotsVarioOut<-paste0(masterDir,"/data_outputs/plots/monthly/variogram")

plotsRFmmOut<-paste0(masterDir,"/data_outputs/plots/monthly/rf_mm")

plotsAnomOut<-paste0(masterDir,"/data_outputs/plots/monthly/anom")

loocvCountyOut<-paste0(masterDir,"/data_outputs/tables/validation/monthly/loocv/county")

metaCountyOut<-paste0(masterDir,"/data_outputs/tables/metadata/monthly/county")

tiffRFmmOut<-paste0(masterDir,"/data_outputs/tiffs/monthly")

statsAllOut<-paste0(masterDir,"/data_outputs/tables/validation/monthly/stats/statewide")

loocvAllOut<-paste0(masterDir,"/data_outputs/tables/validation/monthly/loocv/statewide")

validAllOut<-paste0(masterDir,"/data_outputs/tables/validation/monthly/validate/statewide")

metaStateOut<-paste0(masterDir,"/data_outputs/tables/metadata/monthly/statewide")





#dependencies 

setwd(codeDir)

source("custom_krige_class_fix_functions_Sean.R") #load custom functions



#add in island rasters masks

setwd(dependDir)#set mean rf raster: local pc

bi_mask<-raster("masks/bi_mask.tif")

mn_mask<-raster("masks/mn_mask.tif")

oa_mask<-raster("masks/oa_mask.tif")

ka_mask<-raster("masks/ka_mask.tif")

statewide_mask<-mosaic(mosaic(bi_mask,mn_mask, fun=max),mosaic(oa_mask,ka_mask, fun=max), fun=max)#make statewide raster mask





#Read in file contain Monthly all RF observations load MonYr rf version

setwd(inDir) #set monthly rf data: local pc

dataDate<-Sys.Date()-1 #yesterday

dataYear<-format(dataDate,"%Y")

print(dataYear)
dataYear<-'2022'

RF_MoYr_filename<-paste0("Statewide_Partial_Filled_Monthly_RF_mm_",dataYear,".csv")

RF_MoYr_ALL <- read.csv(file = RF_MoYr_filename , header = TRUE)

head(RF_MoYr_ALL)



#### define monYr col ####

date_cols<-grep("X",names(RF_MoYr_ALL))

date_col_name<-format(dataDate,"X%Y.%m")                
print(date_col_name)
date_col_name<-'X2022.12'
date_col<-grep(date_col_name,names(RF_MoYr_ALL))

s<-Sys.time()#time process



#sub month rf col

data_mon_yr<-gsub("\\.","_",gsub("X","",date_col_name)) #get month year name

RF_MoYr<-RF_MoYr_ALL[!is.na(RF_MoYr_ALL[,date_col]),c(1:(min(date_cols)-1),date_col)]

names(RF_MoYr)[grep("X",names(RF_MoYr))]<-"total_rf_mm" #rename date col as total_rf_mm

#head(RF_MoYr)

  

#define data version

diffDays<-as.numeric(difftime(Sys.Date(),dataDate,units="days"))

if(diffDays==1){version<-"preliminary"

  }else if(diffDays>1 & diffDays<=32){version<-"provisional"

  }else if(diffDays>30 & diffDays<=180){version<-"operational"

    }else if(diffDays>180){version<-"archival"}

  

#Read in Mean RF Grids from a folder and convert to raster 

setwd(paste0(dependDir,"/rf_tiffs"))#set mean rf raster

rf_ras_name<-paste0("staterf_mm",substrRight(data_mon_yr,2),".tif") #make monYr raster name

Mean_RF <- raster(rf_ras_name)#add mean rf of intended month

  

#make spatial dataframe with Lat and Lon

RF_MoYr$x<-RF_MoYr$LON

RF_MoYr$y<-RF_MoYr$LAT

coordinates(RF_MoYr) <- ~x+y

crs(RF_MoYr)<-crs(Mean_RF) #add crs from raster

print(RF_MoYr)

print("Stations Spatial DF Made!")

  

#extract Values from raster

#Creates a df where each column is a daily map with the values extracted for each coordinate

RF_MoYr$RF_Mean_Extract <- extract(Mean_RF, RF_MoYr)

RF_MoYr<-RF_MoYr[!is.na(RF_MoYr$RF_Mean_Extract),] #remove NA stations

#head(RF_MoYr) 

print("Mean Monthly RF extracted!")

  

###Set up for interpolation###

#Calculate anomalies and log anomalies 

C<-1 #define constant C for log + C transformation

RF_MoYr[RF_MoYr$RF_Mean_Extract==0,"RF_Mean_Extract"] <- 0.001 #zeros exsist be sure to not divide by them to avoide dividing into zero during the anomaly calcualtion

RF_MoYr$RF_MoYr_Anom <- RF_MoYr$total_rf_mm / RF_MoYr$RF_Mean_Extract #calc anom" rf obs/mean monthly

RF_MoYr$RF_MoYr_Anom_logK <- log((RF_MoYr$total_rf_mm + C * RF_MoYr$RF_Mean_Extract)/(RF_MoYr$RF_Mean_Extract))#Make suggested transformation a column to back transform

head(RF_MoYr)

print("Anoms and log + C anoms Calculated!")

  

#remove dup/grouped station locations by network priority and dist in meter buffer 'mBuff'

mBuff<-250 #buff dist m

net_priority<-c("HaleNet","CraterNet","Little,HaleNet","HavoNet","HIPPNET","USCRN","ESRL/GMD","NWS","USGS","HydroNet-UaNet","SCAN","NREL","RAWS","unknown","STATE","NREM","PrivateObs","HC&S","COOP","CoCoRaHS","CoCoRaHs")

distbuff<-disaggregate(buffer(RF_MoYr,mBuff))

distbuff$proxID<-as.numeric(row.names(distbuff))

proxID_ras<-rasterize(distbuff, statewide_mask, distbuff$proxID) 

RF_MoYr$proxID<-extract(proxID_ras,RF_MoYr)

RF_MoYr<-RF_MoYr[order(match(RF_MoYr$Network, net_priority)),]

RF_MoYr$map_data<-!duplicated(RF_MoYr$proxID)

print("Close stations removed!")

  

#Write Anom File

setwd(krigInputOut)#export table WD: local PC

rf_month_anom_name<-paste0(data_mon_yr,"_statewide_rf_krig_input.csv") #make anom name

write.csv(RF_MoYr,rf_month_anom_name, row.names = FALSE)

stateSubCounty(as.data.frame(RF_MoYr),rf_month_anom_name,krigInputOutCo) #

print("Anom File Written!")

  

### best krig per county/island(s) ###

##BI best krige

co<-"bi"

bi_kriges_all<-list()

bi_map_check_all<-data.frame()

  

##free all & low fix

bi_varios_low<-varioMaker(anom_sp=RF_MoYr,rasmask=bi_mask,county=co,level="low",none=T)

bi_kriges<-uglyMapFix(anom_sp=RF_MoYr,moYr=data_mon_yr,rasmask=bi_mask,county=co,varioList=bi_varios_low,level="low")

bi_map_check<-bi_kriges[["rf_map_check_low"]]

best_vario_name_bi<-vario_name_func(map_check=bi_map_check)

#check/define best low

if(length(best_vario_name_bi)>0){

	    bi_map_check_all<-rbind(bi_map_check_all,bi_map_check)

    best_krige_list_bi<-bi_kriges[[best_vario_name_bi]]

        print("bi low fixed!")

      } else { 

	          ##med fix

	          bi_map_check_all<-rbind(bi_map_check_all,bi_map_check)

          print("bi med fix check...")

	      bi_varios_med<-varioMaker(anom_sp=RF_MoYr,rasmask=bi_mask,county=co,level="med",none=F)

	      bi_kriges<-append(bi_kriges,uglyMapFix(anom_sp=RF_MoYr,moYr=data_mon_yr,rasmask=bi_mask,county=co,varioList=bi_varios_med,level="med"))

	          bi_map_check<-bi_kriges[["rf_map_check_med"]]

	          best_vario_name_bi<-vario_name_func(map_check=bi_map_check)

		      #check/define best med

		      if(length(best_vario_name_bi)>0){

			            bi_map_check_all<-rbind(bi_map_check_all,bi_map_check)

		        best_krige_list_bi<-bi_kriges[[best_vario_name_bi]]

			      print("bi med fixed!")

			    }else{ 

				          ##high fix

				          bi_map_check_all<-rbind(bi_map_check_all,bi_map_check)

			          print("bi high fix check...")

				        bi_varios_high<-varioMaker(anom_sp=RF_MoYr,rasmask=bi_mask,county=co,level="high",none=F)

				        bi_kriges<-append(bi_kriges,uglyMapFix(anom_sp=RF_MoYr,moYr=data_mon_yr,rasmask=bi_mask,county=co,varioList=bi_varios_high,level="high"))

					      bi_map_check<-bi_kriges[["rf_map_check_high"]]

					      best_vario_name_bi<-vario_name_func(map_check=bi_map_check)

					            #check/define best high

					            if(length(best_vario_name_bi)>0){

							            bi_map_check_all<-rbind(bi_map_check_all,bi_map_check)

					              best_krige_list_bi<-bi_kriges[[best_vario_name_bi]]

						              print("bi high fixed!")

						            } else {

								            #no fix

								            bi_map_check_all<-rbind(bi_map_check_all,bi_map_check)

							            print("bi not fixed... selecting least bad map")

								            best_vario_name_bi<-vario_name_func(bi_map_check_all,noFix=T)

								            best_krige_list_bi<-bi_kriges[[best_vario_name_bi]]

									          }}}#end all fix if else chain

  

##note best vario on map check stats

bi_map_check_all$bestvario<-bi_map_check_all$vario==best_vario_name_bi



##get items from lists 

#get rasters

bi_rf_mm_ras<-best_krige_list_bi$rf_mm_ras

bi_rf_mm_SE_ras<-best_krige_list_bi$rf_mm_SE_ras

bi_krig_anom_ras<-best_krige_list_bi$rf_anom_ras

bi_krig_anom_SE_ras<-best_krige_list_bi$rf_anom_SE_ras

  

##save jpgs

setwd(plotsOut)#export table WD: local PC

  

#multiplot (with function)

jpeg(paste0(plotsMultiOut,"/county/",toupper(co),"/",co,"_",data_mon_yr,"_rfCombo_rf_mm.jpg"),width = 7, height = 7, unit="in", res=300)

multiKrigePlot(kriges=bi_kriges,map_check_all=bi_map_check_all,best_vario_name=best_vario_name_bi)

dev.off()

  

  #variogram

jpeg(paste0(plotsVarioOut,"/county/",toupper(co),"/",co,"_",data_mon_yr,"_vg.jpg"),width = 7, height = 5, unit="in", res=300)

plot(best_krige_list_bi$variogram,sub=paste(toupper(co),data_mon_yr,best_vario_name_bi))

dev.off()

 

#rf plot

jpeg(paste0(plotsRFmmOut,"/county/",toupper(co),"/",co,"_",data_mon_yr,"_rf_mm.jpg"),width = 5, height = 5, unit="in", res=300)

plot(best_krige_list_bi$rf_mm_ras,col=rainbow(300,end=0.8),main=paste(best_vario_name_bi,data_mon_yr,"rainfall mm"), sub=paste0("Class:",best_krige_list_bi$map_validation_df$Class,";  probUgly:",round(best_krige_list_bi$map_validation_df$probUgly*100,1),"%;  Rsq:",round(best_krige_list_bi$map_validation_df$rsq_rf_mm,3),";  RMSE:",round(best_krige_list_bi$map_validation_df$rmse_rf_mm,1)))

dev.off()

  

#anom plot

jpeg(paste0(plotsAnomOut,"/county/",toupper(co),"/",co,"_",data_mon_yr,"_anom.jpg"),width = 5, height = 5, unit="in", res=300)

plot(best_krige_list_bi$rf_anom_ras,col=rev(topo.colors(300)),main=paste(best_vario_name_bi,data_mon_yr,"rainfall anomaly"), sub=paste0("Class:",best_krige_list_bi$map_validation_df$Class,";  probUgly:",round(best_krige_list_bi$map_validation_df$probUgly*100,1),"%;  Rsq:",round(best_krige_list_bi$map_validation_df$rsq_rf_anom,3),";  RMSE:",round(best_krige_list_bi$map_validation_df$rmse_rf_anom,1)))

dev.off()

  

##get loocv, class/error/acc, map stats & station stats tables

bi_loocv_df<-best_krige_list_bi$loocv_df

bi_map_validation_df<-best_krige_list_bi$map_validation_df

  

##get point and raster stats

bi_anom_sp<-best_krige_list_bi$bi_anom_sp

bi_rf_stats<-rbind(rf_dat_stats(dat=bi_anom_sp$total_rf_mm,space=co,var="rf_mm_sta",monyr=data_mon_yr,set=version),

		                      rf_dat_stats(dat=bi_rf_mm_ras,space=co,var="rf_mm_krig",monyr=data_mon_yr,set=version),

				                         rf_dat_stats(dat=bi_rf_mm_SE_ras,space=co,var="rf_mm_SE_krig",monyr=data_mon_yr,set=version))

  

#file names

bi_stations_filename<-gsub("statewide",co,rf_month_anom_name)

bi_rf_ras_filename<-paste0(data_mon_yr,"_",co,"_rf_mm.tif")#rf grid filename

bi_rf_SE_ras_filename<-paste0(data_mon_yr,"_",co,"_rf_mm_SE.tif")#rf SE grid filename

bi_krig_anom_filename<-paste0(data_mon_yr,"_",co,"_anom.tif")

bi_krig_anom_SE_filename<-paste0(data_mon_yr,"_",co,"_anom_SE.tif")

bi_filenames<-c(bi_stations_filename,bi_rf_ras_filename,bi_rf_SE_ras_filename,bi_krig_anom_filename,bi_krig_anom_SE_filename)

  

#write county rf loocv file

setwd(paste0(loocvCountyOut,"/",toupper(co)))#export table WD: local PC

write.csv(bi_loocv_df,paste0(data_mon_yr,"_",co,"_rfmm_loocv.csv"),row.names=F)

print(paste(paste0(data_mon_yr,"_",co,"_rfmm_loocv.csv"),"...written!"))



#write bi meta data file

setwd(paste0(metaCountyOut,"/",toupper(co)))#county meta dir

bi_meta<-metamaker(map_validation_df=bi_map_validation_df,grid=bi_rf_mm_ras,filenames=bi_filenames,datatype=version)

bi_meta_name<-paste0(data_mon_yr,"_",co,"_rf_mm_meta.txt") #rf meta filename

write.table(bi_meta, bi_meta_name,sep ="\t",row.names = F, quote = F)

print(paste(bi_meta_name,"...written!"))





##MN best krige##

co<-"mn"

mn_kriges_all<-list()

mn_map_check_all<-data.frame()



##free all & low fix

mn_varios_low<-varioMaker(anom_sp=RF_MoYr,rasmask=mn_mask,county=co,level="low",none=T)

mn_kriges<-uglyMapFix(anom_sp=RF_MoYr,moYr=data_mon_yr,rasmask=mn_mask,county=co,varioList=mn_varios_low,level="low")

mn_map_check<-mn_kriges[["rf_map_check_low"]]

best_vario_name_mn<-vario_name_func(map_check=mn_map_check)

#check/define best low

if(length(best_vario_name_mn)>0){

	  mn_map_check_all<-rbind(mn_map_check_all,mn_map_check)

  best_krige_list_mn<-mn_kriges[[best_vario_name_mn]]

    print("mn low fixed!")

} else { 

	  ##med fix

	  mn_map_check_all<-rbind(mn_map_check_all,mn_map_check)

  print("mn med fix check...")

    mn_varios_med<-varioMaker(anom_sp=RF_MoYr,rasmask=mn_mask,county=co,level="med",none=F)

    mn_kriges<-append(mn_kriges,uglyMapFix(anom_sp=RF_MoYr,moYr=data_mon_yr,rasmask=mn_mask,county=co,varioList=mn_varios_med,level="med"))

      mn_map_check<-mn_kriges[["rf_map_check_med"]]

      best_vario_name_mn<-vario_name_func(map_check=mn_map_check)

        #check/define best med

        if(length(best_vario_name_mn)>0){

		    mn_map_check_all<-rbind(mn_map_check_all,mn_map_check)

          best_krige_list_mn<-mn_kriges[[best_vario_name_mn]]

	      print("mn med fixed!")

	    }else{ 

		        ##high fix

		        mn_map_check_all<-rbind(mn_map_check_all,mn_map_check)

	        print("mn high fix check...")

		    mn_varios_high<-varioMaker(anom_sp=RF_MoYr,rasmask=mn_mask,county=co,level="high",none=F)

		    mn_kriges<-append(mn_kriges,uglyMapFix(anom_sp=RF_MoYr,moYr=data_mon_yr,rasmask=mn_mask,county=co,varioList=mn_varios_high,level="high"))

		        mn_map_check<-mn_kriges[["rf_map_check_high"]]

		        best_vario_name_mn<-vario_name_func(map_check=mn_map_check)

			    #check/define best high

			    if(length(best_vario_name_mn)>0){

				          mn_map_check_all<-rbind(mn_map_check_all,mn_map_check)

			      best_krige_list_mn<-mn_kriges[[best_vario_name_mn]]

			            print("mn high fixed!")

			          } else {

					        #no fix

					        mn_map_check_all<-rbind(mn_map_check_all,mn_map_check)

				        print("mn not fixed... selecting least bad map")

					      best_vario_name_mn<-vario_name_func(mn_map_check_all,noFix=T)

					      best_krige_list_mn<-mn_kriges[[best_vario_name_mn]]

					          }}}#end all fix if else chain



##note best vario on map check stats

mn_map_check_all$bestvario<-mn_map_check_all$vario==best_vario_name_mn



##get items from lists 

#get rasters

mn_rf_mm_ras<-best_krige_list_mn$rf_mm_ras

mn_rf_mm_SE_ras<-best_krige_list_mn$rf_mm_SE_ras

mn_krig_anom_ras<-best_krige_list_mn$rf_anom_ras

mn_krig_anom_SE_ras<-best_krige_list_mn$rf_anom_SE_ras



##save jpgs

setwd(plotsOut)#export table WD: local PC



#multiplot (with function)

jpeg(paste0(plotsMultiOut,"/county/",toupper(co),"/",co,"_",data_mon_yr,"_rfCombo_rf_mm.jpg"),width = 7, height = 7, unit="in", res=300)

multiKrigePlot(kriges=mn_kriges,map_check_all=mn_map_check_all,best_vario_name=best_vario_name_mn)

dev.off()



#variogram

jpeg(paste0(plotsVarioOut,"/county/",toupper(co),"/",co,"_",data_mon_yr,"_vg.jpg"),width = 7, height = 5, unit="in", res=300)

plot(best_krige_list_mn$variogram,sub=paste(toupper(co),data_mon_yr,best_vario_name_mn))

dev.off()



#rf plot

jpeg(paste0(plotsRFmmOut,"/county/",toupper(co),"/",co,"_",data_mon_yr,"_rf_mm.jpg"),width = 5, height = 5, unit="in", res=300)

plot(best_krige_list_mn$rf_mm_ras,col=rainbow(300,end=0.8),main=paste(best_vario_name_mn,data_mon_yr,"rainfall mm"), sub=paste0("Class:",best_krige_list_mn$map_validation_df$Class,";  probUgly:",round(best_krige_list_mn$map_validation_df$probUgly*100,1),"%;  Rsq:",round(best_krige_list_mn$map_validation_df$rsq_rf_mm,3),";  RMSE:",round(best_krige_list_mn$map_validation_df$rmse_rf_mm,1)))

dev.off()



#anom plot

jpeg(paste0(plotsAnomOut,"/county/",toupper(co),"/",co,"_",data_mon_yr,"_anom.jpg"),width = 5, height = 5, unit="in", res=300)

plot(best_krige_list_mn$rf_anom_ras,col=rev(topo.colors(300)),main=paste(best_vario_name_mn,data_mon_yr,"rainfall anomaly"), sub=paste0("Class:",best_krige_list_mn$map_validation_df$Class,";  probUgly:",round(best_krige_list_mn$map_validation_df$probUgly*100,1),"%;  Rsq:",round(best_krige_list_mn$map_validation_df$rsq_rf_anom,3),";  RMSE:",round(best_krige_list_mn$map_validation_df$rmse_rf_anom,1)))

dev.off()



##get loocv, class/error/acc, map stats & station stats tables

mn_loocv_df<-best_krige_list_mn$loocv_df

mn_map_validation_df<-best_krige_list_mn$map_validation_df



##get point and raster stats

mn_anom_sp<-best_krige_list_mn$mn_anom_sp

mn_rf_stats<-rbind(rf_dat_stats(dat=mn_anom_sp$total_rf_mm,space=co,var="rf_mm_sta",monyr=data_mon_yr,set=version),

		                      rf_dat_stats(dat=mn_rf_mm_ras,space=co,var="rf_mm_krig",monyr=data_mon_yr,set=version),

				                         rf_dat_stats(dat=mn_rf_mm_SE_ras,space=co,var="rf_mm_SE_krig",monyr=data_mon_yr,set=version))



#file names

mn_stations_filename<-gsub("statewide",co,rf_month_anom_name)

mn_rf_ras_filename<-paste0(data_mon_yr,"_",co,"_rf_mm.tif")#rf grid filename

mn_rf_SE_ras_filename<-paste0(data_mon_yr,"_",co,"_rf_mm_SE.tif")#rf SE grid filename

mn_krig_anom_filename<-paste0(data_mon_yr,"_",co,"_anom.tif")

mn_krig_anom_SE_filename<-paste0(data_mon_yr,"_",co,"_anom_SE.tif")

mn_filenames<-c(mn_stations_filename,mn_rf_ras_filename,mn_rf_SE_ras_filename,mn_krig_anom_filename,mn_krig_anom_SE_filename)



#write county rf loocv file

setwd(paste0(loocvCountyOut,"/",toupper(co)))#export table WD: local PC

write.csv(mn_loocv_df,paste0(data_mon_yr,"_",co,"_rfmm_loocv.csv"),row.names=F)

print(paste(paste0(data_mon_yr,"_",co,"_rfmm_loocv.csv"),"...written!"))





#write mn meta data file

setwd(paste0(metaCountyOut,"/",toupper(co)))#county meta dir

mn_meta<-metamaker(map_validation_df=mn_map_validation_df,grid=mn_rf_mm_ras,filenames=mn_filenames,datatype=version)

mn_meta_name<-paste0(data_mon_yr,"_",co,"_rf_mm_meta.txt") #rf meta filename

write.table(mn_meta, mn_meta_name,sep ="\t",row.names = F, quote = F)

print(paste(mn_meta_name,"...written!"))







##OA best krige##

co<-"oa"

oa_kriges_all<-list()

oa_map_check_all<-data.frame()



##free all & low fix

oa_varios_low<-varioMaker(anom_sp=RF_MoYr,rasmask=oa_mask,county=co,level="low",none=T)

oa_kriges<-uglyMapFix(anom_sp=RF_MoYr,moYr=data_mon_yr,rasmask=oa_mask,county=co,varioList=oa_varios_low,level="low")

oa_map_check<-oa_kriges[["rf_map_check_low"]]

best_vario_name_oa<-vario_name_func(map_check=oa_map_check)

#check/define best low

if(length(best_vario_name_oa)>0){

	  oa_map_check_all<-rbind(oa_map_check_all,oa_map_check)

  best_krige_list_oa<-oa_kriges[[best_vario_name_oa]]

    print("oa low fixed!")

} else { 

	  ##med fix

	  oa_map_check_all<-rbind(oa_map_check_all,oa_map_check)

  print("oa med fix check...")

    oa_varios_med<-varioMaker(anom_sp=RF_MoYr,rasmask=oa_mask,county=co,level="med",none=F)

    oa_kriges<-append(oa_kriges,uglyMapFix(anom_sp=RF_MoYr,moYr=data_mon_yr,rasmask=oa_mask,county=co,varioList=oa_varios_med,level="med"))

      oa_map_check<-oa_kriges[["rf_map_check_med"]]

      best_vario_name_oa<-vario_name_func(map_check=oa_map_check)

        #check/define best med

        if(length(best_vario_name_oa)>0){

		    oa_map_check_all<-rbind(oa_map_check_all,oa_map_check)

          best_krige_list_oa<-oa_kriges[[best_vario_name_oa]]

	      print("oa med fixed!")

	    }else{ 

		        ##high fix

		        oa_map_check_all<-rbind(oa_map_check_all,oa_map_check)

	        print("oa high fix check...")

		    oa_varios_high<-varioMaker(anom_sp=RF_MoYr,rasmask=oa_mask,county=co,level="high",none=F)

		    oa_kriges<-append(oa_kriges,uglyMapFix(anom_sp=RF_MoYr,moYr=data_mon_yr,rasmask=oa_mask,county=co,varioList=oa_varios_high,level="high"))

		        oa_map_check<-oa_kriges[["rf_map_check_high"]]

		        best_vario_name_oa<-vario_name_func(map_check=oa_map_check)

			    #check/define best high

			    if(length(best_vario_name_oa)>0){

				          oa_map_check_all<-rbind(oa_map_check_all,oa_map_check)

			      best_krige_list_oa<-oa_kriges[[best_vario_name_oa]]

			            print("oa high fixed!")

			          } else {

					        #no fix

					        oa_map_check_all<-rbind(oa_map_check_all,oa_map_check)

				        print("oa not fixed... selecting least bad map")

					      best_vario_name_oa<-vario_name_func(oa_map_check_all,noFix=T)

					      best_krige_list_oa<-oa_kriges[[best_vario_name_oa]]

					          }}}#end all fix if else chain



##note best vario on map check stats

oa_map_check_all$bestvario<-oa_map_check_all$vario==best_vario_name_oa



##get items from lists 

#get rasters

oa_rf_mm_ras<-best_krige_list_oa$rf_mm_ras

oa_rf_mm_SE_ras<-best_krige_list_oa$rf_mm_SE_ras

oa_krig_anom_ras<-best_krige_list_oa$rf_anom_ras

oa_krig_anom_SE_ras<-best_krige_list_oa$rf_anom_SE_ras



##save jpgs

setwd(plotsOut)#export table WD: local PC



#multiplot (with function)

jpeg(paste0(plotsMultiOut,"/county/",toupper(co),"/",co,"_",data_mon_yr,"_rfCombo_rf_mm.jpg"),width = 7, height = 7, unit="in", res=300)

multiKrigePlot(kriges=oa_kriges,map_check_all=oa_map_check_all,best_vario_name=best_vario_name_oa)

dev.off()



#variogram

jpeg(paste0(plotsVarioOut,"/county/",toupper(co),"/",co,"_",data_mon_yr,"_vg.jpg"),width = 7, height = 5, unit="in", res=300)

plot(best_krige_list_oa$variogram,sub=paste(toupper(co),data_mon_yr,best_vario_name_oa))

dev.off()



#rf plot

jpeg(paste0(plotsRFmmOut,"/county/",toupper(co),"/",co,"_",data_mon_yr,"_rf_mm.jpg"),width = 5, height = 5, unit="in", res=300)

plot(best_krige_list_oa$rf_mm_ras,col=rainbow(300,end=0.8),main=paste(best_vario_name_oa,data_mon_yr,"rainfall mm"), sub=paste0("Class:",best_krige_list_oa$map_validation_df$Class,";  probUgly:",round(best_krige_list_oa$map_validation_df$probUgly*100,1),"%;  Rsq:",round(best_krige_list_oa$map_validation_df$rsq_rf_mm,3),";  RMSE:",round(best_krige_list_oa$map_validation_df$rmse_rf_mm,1)))

dev.off()



#anom plot

jpeg(paste0(plotsAnomOut,"/county/",toupper(co),"/",co,"_",data_mon_yr,"_anom.jpg"),width = 5, height = 5, unit="in", res=300)

plot(best_krige_list_oa$rf_anom_ras,col=rev(topo.colors(300)),main=paste(best_vario_name_oa,data_mon_yr,"rainfall anomaly"), sub=paste0("Class:",best_krige_list_oa$map_validation_df$Class,";  probUgly:",round(best_krige_list_oa$map_validation_df$probUgly*100,1),"%;  Rsq:",round(best_krige_list_oa$map_validation_df$rsq_rf_anom,3),";  RMSE:",round(best_krige_list_oa$map_validation_df$rmse_rf_anom,1)))

dev.off()



##get loocv, class/error/acc, map stats & station stats tables

oa_loocv_df<-best_krige_list_oa$loocv_df

oa_map_validation_df<-best_krige_list_oa$map_validation_df



##get point and raster stats

oa_anom_sp<-best_krige_list_oa$oa_anom_sp

oa_rf_stats<-rbind(rf_dat_stats(dat=oa_anom_sp$total_rf_mm,space=co,var="rf_mm_sta",monyr=data_mon_yr,set=version),

		                      rf_dat_stats(dat=oa_rf_mm_ras,space=co,var="rf_mm_krig",monyr=data_mon_yr,set=version),

				                         rf_dat_stats(dat=oa_rf_mm_SE_ras,space=co,var="rf_mm_SE_krig",monyr=data_mon_yr,set=version))





#file names

oa_stations_filename<-gsub("statewide",co,rf_month_anom_name)

oa_rf_ras_filename<-paste0(data_mon_yr,"_",co,"_rf_mm.tif")#rf grid filename

oa_rf_SE_ras_filename<-paste0(data_mon_yr,"_",co,"_rf_mm_SE.tif")#rf SE grid filename

oa_krig_anom_filename<-paste0(data_mon_yr,"_",co,"_anom.tif")

oa_krig_anom_SE_filename<-paste0(data_mon_yr,"_",co,"_anom_SE.tif")

oa_filenames<-c(oa_stations_filename,oa_rf_ras_filename,oa_rf_SE_ras_filename,oa_krig_anom_filename,oa_krig_anom_SE_filename)



#write county rf loocv file

setwd(paste0(loocvCountyOut,"/",toupper(co)))#export table WD: local PC

write.csv(oa_loocv_df,paste0(data_mon_yr,"_",co,"_rfmm_loocv.csv"),row.names=F)

print(paste(paste0(data_mon_yr,"_",co,"_rfmm_loocv.csv"),"...written!"))





#write oa meta data file

setwd(paste0(metaCountyOut,"/",toupper(co)))#county meta dir

oa_meta<-metamaker(map_validation_df=oa_map_validation_df,grid=oa_rf_mm_ras,filenames=oa_filenames,datatype=version)

oa_meta_name<-paste0(data_mon_yr,"_",co,"_rf_mm_meta.txt") #rf meta filename

write.table(oa_meta, oa_meta_name,sep ="\t",row.names = F, quote = F)

print(paste(oa_meta_name,"...written!"))





##KA best krige##

co<-"ka"

ka_kriges_all<-list()

ka_map_check_all<-data.frame()



##free all & low fix

ka_varios_low<-varioMaker(anom_sp=RF_MoYr,rasmask=ka_mask,county=co,level="low",none=T)

ka_kriges<-uglyMapFix(anom_sp=RF_MoYr,moYr=data_mon_yr,rasmask=ka_mask,county=co,varioList=ka_varios_low,level="low")

ka_map_check<-ka_kriges[["rf_map_check_low"]]

best_vario_name_ka<-vario_name_func(map_check=ka_map_check)

#check/define best low

if(length(best_vario_name_ka)>0){

	  ka_map_check_all<-rbind(ka_map_check_all,ka_map_check)

  best_krige_list_ka<-ka_kriges[[best_vario_name_ka]]

    print("ka low fixed!")

} else { 

	  ##med fix

	  ka_map_check_all<-rbind(ka_map_check_all,ka_map_check)

  print("ka med fix check...")

    ka_varios_med<-varioMaker(anom_sp=RF_MoYr,rasmask=ka_mask,county=co,level="med",none=F)

    ka_kriges<-append(ka_kriges,uglyMapFix(anom_sp=RF_MoYr,moYr=data_mon_yr,rasmask=ka_mask,county=co,varioList=ka_varios_med,level="med"))

      ka_map_check<-ka_kriges[["rf_map_check_med"]]

      best_vario_name_ka<-vario_name_func(map_check=ka_map_check)

        #check/define best med

        if(length(best_vario_name_ka)>0){

		    ka_map_check_all<-rbind(ka_map_check_all,ka_map_check)

          best_krige_list_ka<-ka_kriges[[best_vario_name_ka]]

	      print("ka med fixed!")

	    }else{ 

		        ##high fix

		        ka_map_check_all<-rbind(ka_map_check_all,ka_map_check)

	        print("ka high fix check...")

		    ka_varios_high<-varioMaker(anom_sp=RF_MoYr,rasmask=ka_mask,county=co,level="high",none=F)

		    ka_kriges<-append(ka_kriges,uglyMapFix(anom_sp=RF_MoYr,moYr=data_mon_yr,rasmask=ka_mask,county=co,varioList=ka_varios_high,level="high"))

		        ka_map_check<-ka_kriges[["rf_map_check_high"]]

		        best_vario_name_ka<-vario_name_func(map_check=ka_map_check)

			    #check/define best high

			    if(length(best_vario_name_ka)>0){

				          ka_map_check_all<-rbind(ka_map_check_all,ka_map_check)

			      best_krige_list_ka<-ka_kriges[[best_vario_name_ka]]

			            print("ka high fixed!")

			          } else {

					        #no fix

					        ka_map_check_all<-rbind(ka_map_check_all,ka_map_check)

				        print("ka not fixed... selecting least bad map")

					      best_vario_name_ka<-vario_name_func(ka_map_check_all,noFix=T)

					      best_krige_list_ka<-ka_kriges[[best_vario_name_ka]]

					          }}}#end all fix if else chain



##note best vario on map check stats

ka_map_check_all$bestvario<-ka_map_check_all$vario==best_vario_name_ka



##get items from lists 

#get rasters

ka_rf_mm_ras<-best_krige_list_ka$rf_mm_ras

ka_rf_mm_SE_ras<-best_krige_list_ka$rf_mm_SE_ras

ka_krig_anom_ras<-best_krige_list_ka$rf_anom_ras

ka_krig_anom_SE_ras<-best_krige_list_ka$rf_anom_SE_ras



##save jpgs

setwd(plotsOut)#export table WD: local PC



#multiplot (with function)

jpeg(paste0("multiplot/county/",toupper(co),"/",co,"_",data_mon_yr,"_rfCombo_rf_mm.jpg"),width = 7, height = 7, unit="in", res=300)

multiKrigePlot(kriges=ka_kriges,map_check_all=ka_map_check_all,best_vario_name=best_vario_name_ka)

dev.off()



#variogram

jpeg(paste0(plotsVarioOut,"/county/",toupper(co),"/",co,"_",data_mon_yr,"_vg.jpg"),width = 7, height = 5, unit="in", res=300)

plot(best_krige_list_ka$variogram,sub=paste(toupper(co),data_mon_yr,best_vario_name_ka))

dev.off()



#rf plot

jpeg(paste0(plotsRFmmOut,"/county/",toupper(co),"/",co,"_",data_mon_yr,"_rf_mm.jpg"),width = 5, height = 5, unit="in", res=300)

plot(best_krige_list_ka$rf_mm_ras,col=rainbow(300,end=0.8),main=paste(best_vario_name_ka,data_mon_yr,"rainfall mm"), sub=paste0("Class:",best_krige_list_ka$map_validation_df$Class,";  probUgly:",round(best_krige_list_ka$map_validation_df$probUgly*100,1),"%;  Rsq:",round(best_krige_list_ka$map_validation_df$rsq_rf_mm,3),";  RMSE:",round(best_krige_list_ka$map_validation_df$rmse_rf_mm,1)))

dev.off()



#anom plot

jpeg(paste0(plotsAnomOut,"/county/",toupper(co),"/",co,"_",data_mon_yr,"_anom.jpg"),width = 5, height = 5, unit="in", res=300)

plot(best_krige_list_ka$rf_anom_ras,col=rev(topo.colors(300)),main=paste(best_vario_name_ka,data_mon_yr,"rainfall anomaly"), sub=paste0("Class:",best_krige_list_ka$map_validation_df$Class,";  probUgly:",round(best_krige_list_ka$map_validation_df$probUgly*100,1),"%;  Rsq:",round(best_krige_list_ka$map_validation_df$rsq_rf_anom,3),";  RMSE:",round(best_krige_list_ka$map_validation_df$rmse_rf_anom,1)))

dev.off()



##get loocv, class/error/acc, map stats & station stats tables

ka_loocv_df<-best_krige_list_ka$loocv_df

ka_map_validation_df<-best_krige_list_ka$map_validation_df



##get point and raster stats

ka_anom_sp<-best_krige_list_ka$ka_anom_sp

ka_rf_stats<-rbind(rf_dat_stats(dat=ka_anom_sp$total_rf_mm,space=co,var="rf_mm_sta",monyr=data_mon_yr,set=version),

		                      rf_dat_stats(dat=ka_rf_mm_ras,space=co,var="rf_mm_krig",monyr=data_mon_yr,set=version),

				                         rf_dat_stats(dat=ka_rf_mm_SE_ras,space=co,var="rf_mm_SE_krig",monyr=data_mon_yr,set=version))





#file names

ka_stations_filename<-gsub("statewide",co,rf_month_anom_name)

ka_rf_ras_filename<-paste0(data_mon_yr,"_",co,"_rf_mm.tif")#rf grid filename

ka_rf_SE_ras_filename<-paste0(data_mon_yr,"_",co,"_rf_mm_SE.tif")#rf SE grid filename

ka_krig_anom_filename<-paste0(data_mon_yr,"_",co,"_anom.tif")

ka_krig_anom_SE_filename<-paste0(data_mon_yr,"_",co,"_anom_SE.tif")

ka_filenames<-c(ka_stations_filename,ka_rf_ras_filename,ka_rf_SE_ras_filename,ka_krig_anom_filename,ka_krig_anom_SE_filename)



#write county rf loocv file

setwd(paste0(loocvCountyOut,"/",toupper(co)))#export table WD: local PC

write.csv(ka_loocv_df,paste0(data_mon_yr,"_",co,"_rfmm_loocv.csv"),row.names=F)

print(paste(paste0(data_mon_yr,"_",co,"_rfmm_loocv.csv"),"...written!"))





#write ka meta data file

setwd(paste0(metaCountyOut,"/",toupper(co)))#county meta dir

ka_meta<-metamaker(map_validation_df=ka_map_validation_df,grid=ka_rf_mm_ras,filenames=ka_filenames,datatype=version)

ka_meta_name<-paste0(data_mon_yr,"_",co,"_rf_mm_meta.txt") #rf meta filename

write.table(ka_meta, ka_meta_name,sep ="\t",row.names = F, quote = F)

print(paste(ka_meta_name,"...written!"))

### County KRIGING maps pau ###  

  

###Final all county plot

setwd(paste0(plotsRFmmOut,"/statewide"))

jpeg(paste0("all_co_",data_mon_yr,"_rf_mm.jpg"),width = 8, height = 6, unit="in", res=300)

par(mfrow=c(2,2))

par(mai = c(0.35,0.35,0.35,0.35))

  plot(bi_rf_mm_ras,col=rainbow(300,end=0.8),axes=FALSE, box=FALSE,horizontal = TRUE,

              legend.width=1,legend.args=list(text='mm', side=2, font=3,cex=0.8))

  title(sub=best_vario_name_bi,line=-0.5)

    title(main=paste("BI:",data_mon_yr),line=-0.3)

    plot(mn_rf_mm_ras,col=rainbow(300,end=0.8),axes=FALSE, box=FALSE,horizontal = TRUE,

	        legend.width=1,legend.args=list(text='mm', side=2, font=3,cex=0.8))

      title(sub=best_vario_name_mn,line=-0.5)

      title(main=paste("MN:",data_mon_yr),line=-0.3)

        plot(oa_rf_mm_ras,col=rainbow(300,end=0.8),

	            axes=FALSE, box=FALSE,horizontal = TRUE,

		           legend.width=1,legend.args=list(text='mm', side=2, font=3,cex=0.8))

        title(sub=best_vario_name_oa,line=-1)

	  title(main=paste("OA:",data_mon_yr),line=-0.3)

	  plot(ka_rf_mm_ras,col=rainbow(300,end=0.8),axes=FALSE, box=FALSE,horizontal = TRUE,

	              legend.width=1,legend.args=list(text='mm', side=2, font=3,cex=0.8))

	    title(sub=best_vario_name_ka,line=-1)

	    title(main=paste("KA:",data_mon_yr),line=-0.3)

	    dev.off()

	      

	    ###write rasters

	    ## write island tif files ### 

	    #RF mm

	    setwd(paste0(tiffRFmmOut,"/rf_mm/county")) #output tiffs dir

	    writeRaster(bi_rf_mm_ras,paste0("BI/",bi_rf_ras_filename),overwrite=TRUE)

	    writeRaster(mn_rf_mm_ras,paste0("MN/",mn_rf_ras_filename),overwrite=TRUE)

	    writeRaster(oa_rf_mm_ras,paste0("OA/",oa_rf_ras_filename),overwrite=TRUE)

	    writeRaster(ka_rf_mm_ras,paste0("KA/",ka_rf_ras_filename),overwrite=TRUE)

	    print("county RF tiffs Written!")

	      

	    #RF SE

	    setwd(paste0(tiffRFmmOut,"/rf_mm_se/county")) #output tiffs dir

	    writeRaster(bi_rf_mm_SE_ras,paste0("BI/",bi_rf_SE_ras_filename),overwrite=TRUE)

	    writeRaster(mn_rf_mm_SE_ras,paste0("MN/",mn_rf_SE_ras_filename),overwrite=TRUE)

	    writeRaster(oa_rf_mm_SE_ras,paste0("OA/",oa_rf_SE_ras_filename),overwrite=TRUE)

	    writeRaster(ka_rf_mm_SE_ras,paste0("KA/",ka_rf_SE_ras_filename),overwrite=TRUE)

	    print("county SE tiffs Written!")

	      

	    #Anoms

	    setwd(paste0(tiffRFmmOut,"/anom/county")) #output tiffs dir

	    writeRaster(bi_krig_anom_ras,paste0("BI/",bi_krig_anom_filename),overwrite=TRUE)

	    writeRaster(mn_krig_anom_ras,paste0("MN/",mn_krig_anom_filename),overwrite=TRUE)

	    writeRaster(oa_krig_anom_ras,paste0("OA/",oa_krig_anom_filename),overwrite=TRUE)

	    writeRaster(ka_krig_anom_ras,paste0("KA/",ka_krig_anom_filename),overwrite=TRUE)

	    print("county Anom tiffs Written!")

	      

	    #Anom SE 

	    setwd(paste0(tiffRFmmOut,"/anom_se/county")) #output tiffs dir

	    writeRaster(bi_krig_anom_SE_ras,paste0("BI/",bi_krig_anom_SE_filename),overwrite=TRUE)

	    writeRaster(mn_krig_anom_SE_ras,paste0("MN/",mn_krig_anom_SE_filename),overwrite=TRUE)

	    writeRaster(oa_krig_anom_SE_ras,paste0("OA/",oa_krig_anom_SE_filename),overwrite=TRUE)

	    writeRaster(ka_krig_anom_SE_ras,paste0("KA/",ka_krig_anom_SE_filename),overwrite=TRUE)

	    print("county SE anom tiffs Written!")

	      

	    #statewide file names

	    hi_rf_ras_filename<-paste0(data_mon_yr,"_statewide_rf_mm.tif")#rf grid filename

	    hi_rf_SE_ras_filename<-paste0(data_mon_yr,"_statewide_rf_mm_SE.tif")#rf SE grid filename

	    hi_krig_anom_filename<-paste0(data_mon_yr,"_statewide_anom.tif")

	    hi_krig_anom_SE_filename<-paste0(data_mon_yr,"_statewide_anom_SE.tif")

	    hi_filenames<-c(rf_month_anom_name,hi_rf_ras_filename,hi_rf_SE_ras_filename,hi_krig_anom_filename,hi_krig_anom_SE_filename)

	      

	    #combind and write statewide data

	    #RF mm

	    setwd(paste0(tiffRFmmOut,"/rf_mm/statewide")) #output tiffs dir

	    hi_statewide_rf_mm_ras<-mosaic(bi_rf_mm_ras, mn_rf_mm_ras, fun=max)

	    hi_statewide_rf_mm_ras<-mosaic(hi_statewide_rf_mm_ras, oa_rf_mm_ras, fun=max)

	    hi_statewide_rf_mm_ras<-mosaic(hi_statewide_rf_mm_ras, ka_rf_mm_ras, fun=max)

	    #plot(hi_statewide_rf_mm_ras)

	    writeRaster(hi_statewide_rf_mm_ras,hi_rf_ras_filename,overwrite=TRUE)

	    print("Statewide RF tiff Written!")

	      

	    #RF mm SE

	    setwd(paste0(tiffRFmmOut,"/rf_mm_se/statewide")) #output tiffs dir

	    hi_statewide_rf_mm_SE_ras<-mosaic(bi_rf_mm_SE_ras, mn_rf_mm_SE_ras, fun=max)

	    hi_statewide_rf_mm_SE_ras<-mosaic(hi_statewide_rf_mm_SE_ras, oa_rf_mm_SE_ras, fun=max)

	    hi_statewide_rf_mm_SE_ras<-mosaic(hi_statewide_rf_mm_SE_ras, ka_rf_mm_SE_ras, fun=max)

	    #plot(hi_statewide_rf_mm_SE_ras)

	    writeRaster(hi_statewide_rf_mm_SE_ras,hi_rf_SE_ras_filename,overwrite=TRUE)

	    print("Statewide SE tiff Written!")

	      

	    #Anom

	    setwd(paste0(tiffRFmmOut,"/anom/statewide")) #output tiffs dir

	    hi_statewide_krig_anom_ras<-mosaic(bi_krig_anom_ras, mn_krig_anom_ras, fun=max)

	    hi_statewide_krig_anom_ras<-mosaic(hi_statewide_krig_anom_ras, oa_krig_anom_ras, fun=max)

	    hi_statewide_krig_anom_ras<-mosaic(hi_statewide_krig_anom_ras, ka_krig_anom_ras, fun=max)

	    #plot(hi_statewide_krig_anom_ras)

	    writeRaster(hi_statewide_krig_anom_ras,hi_krig_anom_filename,overwrite=TRUE)

	    print("Statewide Anom tiff Written!")

	      

	    #Anom SE

	    setwd(paste0(tiffRFmmOut,"/anom_se/statewide")) #output tiffs dir

	    hi_statewide_anom_SE_ras<-mosaic(bi_krig_anom_SE_ras, mn_krig_anom_SE_ras, fun=max)

	    hi_statewide_anom_SE_ras<-mosaic(hi_statewide_anom_SE_ras, oa_krig_anom_SE_ras, fun=max)

	    hi_statewide_anom_SE_ras<-mosaic(hi_statewide_anom_SE_ras, ka_krig_anom_SE_ras, fun=max)

	    #plot(hi_statewide_anom_SE_ras)

	    writeRaster(hi_statewide_anom_SE_ras,hi_krig_anom_SE_filename,overwrite=TRUE)

	    print("Statewide Anom SE tiff Written!")

	    #end write rasters

	      

	    ###combind and write tables and final statewide meta

	    #rf stats

	    setwd(statsAllOut)#export table WD

	    all_rf_stats<-rbind(bi_rf_stats,mn_rf_stats,oa_rf_stats,ka_rf_stats)

	    write.csv(all_rf_stats,paste0(data_mon_yr,"_all_rf_stats.csv"),row.names=F)

	      

	    #rf loocv

	    setwd(loocvAllOut)#export table WD: local PC

	    all_loocv_df<-rbind(bi_loocv_df,mn_loocv_df,oa_loocv_df,ka_loocv_df)

	    write.csv(all_loocv_df,paste0(data_mon_yr,"_all_rfmm_loocv.csv"),row.names=F)

	      

	    #rf validations

	    setwd(validAllOut)#export table WD: local PC

	    all_map_validation_df<-rbind(bi_map_check_all,mn_map_check_all,oa_map_check_all,ka_map_check_all)

	    write.csv(all_map_validation_df,paste0("all/",data_mon_yr,"_all_map_validation.csv"),row.names=F)

	    best_map_validation_df<-rbind(bi_map_validation_df,mn_map_validation_df,oa_map_validation_df,ka_map_validation_df)

	    write.csv(best_map_validation_df,paste0("best/",data_mon_yr,"_best_map_validation.csv"),row.names=F)



	    #write final statewide metadata

	    setwd(metaStateOut)

	    hi_meta<-metamaker(map_validation_df=best_map_validation_df,grid=hi_statewide_rf_mm_ras,filenames=hi_filenames,datatype=version,statewide=T,loocv_df=all_loocv_df)

	    hi_meta_name<-paste0(data_mon_yr,"_statewide_rf_mm_meta.txt")#rf meta filename

	    write.table(hi_meta, hi_meta_name ,sep ="\t",row.names = F, quote = F)



	    e<-Sys.time()

	    print(e-s)# time taken



	    #code pau
