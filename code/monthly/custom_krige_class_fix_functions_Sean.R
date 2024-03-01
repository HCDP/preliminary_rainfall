

#packages used

#install.packages(c("dplyr","randomForest","raster","automap","gstat","Metrics","svMisc"))





#custom funtions to for auto kriging 



#substrRight func 

substrRight <- function(x, n){

	  substr(x, nchar(x)-n+1, nchar(x))

  }#end substrRight function



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



#rf_dat_stats function

rf_dat_stats<-function(dat,space,var,monyr,set){ #fun to calc stats per dataset

		if(!is.vector(dat)){

				 vals<-as.numeric(na.omit(values(dat)))

	 stations<-as.numeric(NA)

	 	 }else{

			 	 vals<-as.numeric(dat)

	 	 stations<-length(dat)

		 	 }

	stats_df<- data.frame(

			      		     var=var,county=toupper(space),

					     		     monYr=monyr,dataset=set,

					     		     min=min(vals),mean=mean(vals),

							     		     median=median(vals),max=max(vals),

							     		     SD=sd(vals),SE=sd(vals)/sqrt(length(vals)),

									     	 	     stations=as.numeric(NA))



		if(!is.na(stations)){stats_df$stations<-stations}

		return(stats_df)

			} #rf_dat_stats func end











#is ugly class/prob function

isUglyClass<-function(rf_points,rf_ras,countyIntials){

	    options(warn = -1)

  #packages

	require(raster)

		require(dplyr)

		require(randomForest)



		  #load rf model

			setwd("/home/hawaii_climate_products_container/preliminary/rainfall/dependencies/monthly")

			uglyClassfier<- readRDS("./isUgly_RF_classfier.rds")



			  #make dist to station rasters

				crs(rf_points)<-crs(rf_ras)

				dist<-mask(distanceFromPoints(rf_ras,rf_points),rf_ras)

					#plot(dist)

					thresh<-dist<1000 #mask < 1000m



				 #make terrain vars

					slope<-terrain(rf_ras, opt="slope", unit="degrees",neigbors=8)

						rough<-terrain(rf_ras, opt="roughness")

						tri<-terrain(rf_ras, opt="tri")

							curve<-terrain(terrain(rf_ras, opt="slope", unit="degrees",neigbors=8), opt="slope", unit="degrees",neigbors=8)



						 #stack layers make df

							ras_stack<-stack(rf_ras,dist,thresh,slope,rough,tri,curve)

								names(ras_stack)<-c("rf_mm","distm","thresh","slopeDeg","rough","tri","curvature")

							 	df_stack<-as.data.frame(ras_stack,xy=T,na.rm=T)



								 #make monYrCo df

									rfdist_df<-data.frame(

											       				county=c(countyIntials,"bi","ka","mn","oa"),

																			rf_neg=c(min(df_stack$rf_mm)<0,"FALSE", "TRUE","FALSE", "TRUE"),

																			rf_min=min(df_stack$rf_mm),

																							rf_max=max(df_stack$rf_mm),

																							rf_mean=mean(df_stack$rf_mm),

																											rf_median=median(df_stack$rf_mm),

																											rf_sd=sd(df_stack$rf_mm),

																															rf_var=var(df_stack$rf_mm),

																															slopeDeg_min=min(df_stack$slopeDeg),

																																			slopeDeg_max=max(df_stack$slopeDeg),

																																			slopeDeg_mean=mean(df_stack$slopeDeg),

																																							slopeDeg_median=median(df_stack$slopeDeg),

																																							slopeDeg_sd=sd(df_stack$slopeDeg),

																																											slopeDeg_var=var(df_stack$slopeDeg),

																																											rough_min=min(df_stack$rough),

																																															rough_max=max(df_stack$rough),

																																															rough_mean=mean(df_stack$rough),

																																																			rough_median=median(df_stack$rough),

																																																			rough_sd=sd(df_stack$rough),

																																																							rough_var=var(df_stack$rough),

																																																							tri_min=min(df_stack$tri),

																																																											tri_max=max(df_stack$tri),

																																																											tri_mean=mean(df_stack$tri),

																																																															tri_median=median(df_stack$tri),

																																																															tri_sd=sd(df_stack$tri),

																																																																			tri_var=var(df_stack$tri),

																																																																			curve_min=min(df_stack$curvature),

																																																																							curve_max=max(df_stack$curvature),

																																																																							curve_mean=mean(df_stack$curvature),

																																																																											curve_median=median(df_stack$curvature),

																																																																											curve_sd=sd(df_stack$curvature),

																																																																															curve_var=var(df_stack$curvature),

																																																																															cor_slp=cor(df_stack$slopeDeg,df_stack$distm),

																																																																																			cor_rough=cor(df_stack$rough,df_stack$distm),

																																																																																			cor_tri=cor(df_stack$tri,df_stack$distm),

																																																																																							cor_curve=cor(df_stack$curvature,df_stack$distm),

																																																																																							rsq_slp=as.numeric(summary(lm(df_stack$slopeDeg~df_stack$distm))[9]),

																																																																																											rsq_rough=as.numeric(summary(lm(df_stack$rough~df_stack$distm))[9]),

																																																																																											rsq_tri=as.numeric(summary(lm(df_stack$tri~df_stack$distm))[9]),

																																																																																															rsq_curve=as.numeric(summary(lm(df_stack$curvature~df_stack$distm))[9]),

																																																																																															pval_slp=c((summary(lm(df_stack$slopeDeg~df_stack$distm))$coefficients[,"Pr(>|t|)"][[2]]<0.05),"FALSE", "TRUE","FALSE", "TRUE"),

																																																																																																			pval_rough=c((summary(lm(df_stack$rough~df_stack$distm))$coefficients[,"Pr(>|t|)"][[2]]<0.05),"FALSE", "TRUE","FALSE", "TRUE"),

																																																																																																			pval_tri=c((summary(lm(df_stack$tri~df_stack$distm))$coefficients[,"Pr(>|t|)"][[2]]<0.05),"FALSE", "TRUE","FALSE", "TRUE"),

																																																																																																							pval_curve=c((summary(lm(df_stack$curvature~df_stack$distm))$coefficients[,"Pr(>|t|)"][[2]]<0.05),"FALSE", "TRUE","FALSE", "TRUE"),

																																																																																																							rf_min_thresh=min(df_stack[as.logical(df_stack$thresh),"rf_mm"]),

																																																																																																											rf_max_thresh=max(df_stack[as.logical(df_stack$thresh),"rf_mm"]),

																																																																																																											rf_mean_thresh=mean(df_stack[as.logical(df_stack$thresh),"rf_mm"]),

																																																																																																															rf_median_thresh=median(df_stack[as.logical(df_stack$thresh),"rf_mm"]),

																																																																																																															rf_sd_thresh=sd(df_stack[as.logical(df_stack$thresh),"rf_mm"]),

																																																																																																																			rf_var_thresh=var(df_stack[as.logical(df_stack$thresh),"rf_mm"]),

																																																																																																																			slopeDeg_min_thresh=min(df_stack[as.logical(df_stack$thresh),"slopeDeg"]),

																																																																																																																							slopeDeg_max_thresh=max(df_stack[as.logical(df_stack$thresh),"slopeDeg"]),

																																																																																																																							slopeDeg_mean_thresh=mean(df_stack[as.logical(df_stack$thresh),"slopeDeg"]),

																																																																																																																											slopeDeg_median_thresh=median(df_stack[as.logical(df_stack$thresh),"slopeDeg"]),

																																																																																																																											slopeDeg_sd_thresh=sd(df_stack[as.logical(df_stack$thresh),"slopeDeg"]),

																																																																																																																															slopeDeg_var_thresh=var(df_stack[as.logical(df_stack$thresh),"slopeDeg"]),

																																																																																																																															rough_min_thresh=min(df_stack[as.logical(df_stack$thresh),"rough"]),

																																																																																																																																			rough_max_thresh=max(df_stack[as.logical(df_stack$thresh),"rough"]),

																																																																																																																																			rough_mean_thresh=mean(df_stack[as.logical(df_stack$thresh),"rough"]),

																																																																																																																																							rough_median_thresh=median(df_stack[as.logical(df_stack$thresh),"rough"]),

																																																																																																																																							rough_sd_thresh=sd(df_stack[as.logical(df_stack$thresh),"rough"]),

																																																																																																																																											rough_var_thresh=var(df_stack[as.logical(df_stack$thresh),"rough"]),

																																																																																																																																											tri_min_thresh=min(df_stack[as.logical(df_stack$thresh),"tri"]),

																																																																																																																																															tri_max_thresh=max(df_stack[as.logical(df_stack$thresh),"tri"]),

																																																																																																																																															tri_mean_thresh=mean(df_stack[as.logical(df_stack$thresh),"tri"]),

																																																																																																																																																			tri_median_thresh=median(df_stack[as.logical(df_stack$thresh),"tri"]),

																																																																																																																																																			tri_sd_thresh=sd(df_stack[as.logical(df_stack$thresh),"tri"]),

																																																																																																																																																							tri_var_thresh=var(df_stack[as.logical(df_stack$thresh),"tri"]),

																																																																																																																																																							curve_min_thresh=min(df_stack[as.logical(df_stack$thresh),"curvature"]),

																																																																																																																																																											curve_max_thresh=max(df_stack[as.logical(df_stack$thresh),"curvature"]),

																																																																																																																																																											curve_mean_thresh=mean(df_stack[as.logical(df_stack$thresh),"curvature"]),

																																																																																																																																																															curve_median_thresh=median(df_stack[as.logical(df_stack$thresh),"curvature"]),

																																																																																																																																																															curve_sd_thresh=sd(df_stack[as.logical(df_stack$thresh),"curvature"]),

																																																																																																																																																																			curve_var_thresh=var(df_stack[as.logical(df_stack$thresh),"curvature"]),

																																																																																																																																																																			cor_slp_thresh=cor(df_stack[as.logical(df_stack$thresh),"slopeDeg"],df_stack[as.logical(df_stack$thresh),"distm"]),

																																																																																																																																																																							cor_rough_thresh=cor(df_stack[as.logical(df_stack$thresh),"rough"],df_stack[as.logical(df_stack$thresh),"distm"]),

																																																																																																																																																																							cor_tri_thresh=cor(df_stack[as.logical(df_stack$thresh),"tri"],df_stack[as.logical(df_stack$thresh),"distm"]),

																																																																																																																																																																											cor_curve_thresh=cor(df_stack[as.logical(df_stack$thresh),"curvature"],df_stack[as.logical(df_stack$thresh),"distm"]),

																																																																																																																																																																											rsq_slp_thresh=as.numeric(summary(lm(df_stack[as.logical(df_stack$thresh),"slopeDeg"]~df_stack[as.logical(df_stack$thresh),"distm"]))[9]),

																																																																																																																																																																															rsq_rough_thresh=as.numeric(summary(lm(df_stack[as.logical(df_stack$thresh),"rough"]~df_stack[as.logical(df_stack$thresh),"distm"]))[9]),

																																																																																																																																																																															rsq_tri_thresh=as.numeric(summary(lm(df_stack[as.logical(df_stack$thresh),"tri"]~df_stack[as.logical(df_stack$thresh),"distm"]))[9]),

																																																																																																																																																																																			rsq_curve_thresh=as.numeric(summary(lm(df_stack[as.logical(df_stack$thresh),"curvature"]~df_stack[as.logical(df_stack$thresh),"distm"]))[9]),

																																																																																																																																																																																			pval_slp_thresh=c((summary(lm(df_stack[as.logical(df_stack$thresh),"slopeDeg"]~df_stack[as.logical(df_stack$thresh),"distm"]))$coefficients[,"Pr(>|t|)"][[2]]<0.05),"FALSE", "TRUE","FALSE", "TRUE"),

																																																																																																																																																																																							pval_rough_thresh=c((summary(lm(df_stack[as.logical(df_stack$thresh),"rough"]~df_stack[as.logical(df_stack$thresh),"distm"]))$coefficients[,"Pr(>|t|)"][[2]]<0.05),"FALSE", "TRUE","FALSE", "TRUE"),

																																																																																																																																																																																							pval_tri_thresh=c((summary(lm(df_stack[as.logical(df_stack$thresh),"tri"]~df_stack[as.logical(df_stack$thresh),"distm"]))$coefficients[,"Pr(>|t|)"][[2]]<0.05),"FALSE", "TRUE","FALSE", "TRUE"),

																																																																																																																																																																																											pval_curve_thresh=c((summary(lm(df_stack[as.logical(df_stack$thresh),"curvature"]~df_stack[as.logical(df_stack$thresh),"distm"]))$coefficients[,"Pr(>|t|)"][[2]]<0.05),"FALSE", "TRUE","FALSE", "TRUE")

																																																																																																																																																																																											)

									#prep df

									rfdist_df <- rfdist_df %>% mutate_if(is.character,as.factor)

										rfdist_df <- rfdist_df %>% mutate_if(is.logical,as.factor)

										rfdist_df <- rfdist_df[1,]

											#return(head(rfdist_df))

											#predict using trained classfier

											isUgly_prediction <- predict(uglyClassfier, newdata=rfdist_df)

											isUgly_prob <- as.numeric(as.data.frame(predict(uglyClassfier, newdata=rfdist_df, type="prob"))[,"isUgly"])

												return(data.frame(Class=as.character(isUgly_prediction),probUgly=as.numeric(isUgly_prob)))



}#isUgly funtion end



varioMaker<-function(anom_sp,rasmask,county,level,none){



			#required packages

			require(automap)

		require(gstat)



				#define values (nug, range & sill) if map is uglyfor each county

				NugRanSil<-read.csv("/home/hawaii_climate_products_container/preliminary/rainfall/dependencies/monthly/krige_county_vars_final.csv",stringsAsFactors=F)

				#print(NugRanSil)



				#subset rf points and make temp points

				anom_sp<-anom_sp[!is.na(extract(rasmask,anom_sp)) & anom_sp$map_data,]#subset county and unique prox

						crs(anom_sp)<-CRS() #remove crs

						temppoints<-SpatialPoints(as.data.frame(rasmask,xy=T,na.rm=T)[,c(1,2)]) #make spatial points data frame with x y coords from mask and remove NA pixels

								#crs(rasmask)<-CRS() #remove crs

								

								#make all vatiograms for level

								varioList<-list()

								if(as.logical(none)){

											 varioList[["vari_free_all"]]<-autofitVariogram(as.formula("RF_MoYr_Anom_logK ~ 1"), anom_sp, model="Ste", )

										}

								 		tryCatch({varioList[[paste0("vari_fix_nug_",level)]]<-autofitVariogram(as.formula("RF_MoYr_Anom_logK ~ 1"), anom_sp, model="Ste", fix.values = as.numeric(c(as.numeric(NugRanSil[NugRanSil$county==toupper(county),c(paste0(level,"_nug"))]),"NA","NA")))

													}, error=function(e){})

										

										tryCatch({varioList[[paste0("vari_fix_rng_",level)]]<-autofitVariogram(as.formula("RF_MoYr_Anom_logK ~ 1"), anom_sp, model="Ste", fix.values = as.numeric(c("NA",as.numeric(NugRanSil[NugRanSil$county==toupper(county),paste0(level,"_rng")]),"NA")))

													}, error=function(e){})



												tryCatch({varioList[[paste0("vari_free_sil_",level)]]<-autofitVariogram(as.formula("RF_MoYr_Anom_logK ~ 1"), anom_sp, model="Ste", fix.values = as.numeric(c(as.numeric(NugRanSil[NugRanSil$county==toupper(county),c(paste0(level,"_nug"),paste0(level,"_rng"))]),"NA")))

													  		}, error=function(e){})



												return(varioList)

														}#end varioMaker function



#ugly map fix function

uglyMapFix<-function(anom_sp,moYr,rasmask,county,varioList,level){



			#required packages

			require(automap)

		require(gstat)

				require(Metrics)

				require(svMisc)



						#define values (nug, range & sill) if map is uglyfor each county

						setwd("/home/hawaii_climate_products_container/preliminary/rainfall/dependencies/monthly")

						NugRanSil<-read.csv("krige_county_vars_final.csv",stringsAsFactors=F)



								#make masterlist to store all the goodies

								runList<-list()

								rf_fix_df_all<-data.frame()#blank df to storm final all vario df



										#subset rf points and make temp points

										anom_sp<-anom_sp[!is.na(extract(rasmask,anom_sp)) & anom_sp$map_data,]#subset county and unique prox

										crs(anom_sp)<-CRS() #remove crs

												temppoints<-SpatialPoints(as.data.frame(rasmask,xy=T,na.rm=T)[,c(1,2)]) #make spatial points data frame with x y coords from mask and remove NA pixels

												

												for(v in 1:length(varioList)){

																vario<-varioList[[v]]

															krigeObj<-krige(RF_MoYr_Anom_logK ~ 1 ,anom_sp, temppoints, model=vario$var_model)

																		krig_logC_anom<-rasterize(krigeObj, rasmask, krigeObj$var1.pred) #make krig log points into raster

																		krig_logC_anom_SE<-rasterize(krigeObj, rasmask, krigeObj$var1.var) #make krig SD logk points into raster

																			

																					#convert to rf mm units

																					rf_mm_ras<- round((exp(krig_logC_anom) - C)*Mean_RF,10) #back transform to get rf mm total: note not Anom

																					rf_mm_SE_ras<- round((exp(krig_logC_anom_SE) - C)*Mean_RF,10) #back transform SE to get rf mm total: note not Anom

																								rf_anom_ras<- rf_mm_ras/(Mean_RF+0.00000001) #back transform to get rf mm total: note not Anom

																								rf_anom_SE_ras<- rf_mm_SE_ras/(Mean_RF+0.00000001) #back transform to get rf mm total: note not Anom

																											

																											#class ugly with prob

																											uglyClass<-isUglyClass(rf_points=anom_sp,rf_ras=rf_mm_ras,countyIntials=county)



																											#save IF neg RF values produced and edit neg RF to zero

																											negRFzeroGrid<-minValue(rf_mm_ras)<0 #save neg value were present

																														krigeRFmin<-minValue(rf_mm_ras)

																														krigeRFmax<-maxValue(rf_mm_ras)

																																	negRFpixCount<-sum(na.omit(values(rf_mm_ras<0)))

																																	gridPixCount<-length(na.omit(values(rf_mm_ras)))

																																				rf_mm_ras[rf_mm_ras<0]<-0 #edit neg values to zero			



																																				#loocv loop

																																				krige_rf<-as.numeric()

																																							sink("deleteMe.txt")#silence and save everything (delete later)

																																							for(k in 1:nrow(anom_sp)){

																																												krige_loocv<-krige(RF_MoYr_Anom_logK ~ 1 ,anom_sp[-k,],anom_sp[k,], model=vario$var_model)

																																											krige_rf<-c(krige_rf,as.numeric(krige_loocv$var1.pred))

																																															}#end loocv

																																							 			sink()#stop message sink

																																										unlink("deleteMe.txt")#delete message file



																																													#get and add rf vals

																																													loocv_df<-as.data.frame(anom_sp)[,c("SKN","Island","Network","total_rf_mm","RF_Mean_Extract","RF_MoYr_Anom")] #make rf sub df

																																													names(loocv_df)<-c("SKN","Island","Network","obs_rf","RF_Mean_Extract","obs_anom")

																																																loocv_df$pred_rf<-round(exp(krige_rf - C)*loocv_df$RF_Mean_Extract,10)#back transform and add RF col

																																																loocv_df$pred_anom<-loocv_df$pred_rf/(loocv_df$RF_Mean_Extract+0.00000001) #calc rf anom

																																																			loocv_df$negRFzero<-loocv_df$pred_rf<0

																																																			loocv_df[loocv_df$negRFzero,"pred_rf"]<-as.numeric(0) #change neg RF pred to zero

																																																						loocv_df$county<-county

																																																						loocv_df$level<-level

																																																									loocv_df$vario<-names(varioList)[v]

																																																									loocv_df<-loocv_df[,c("SKN","Island","county","Network","vario","RF_Mean_Extract","level","obs_rf","obs_anom","pred_rf","pred_anom","negRFzero")]

																																																												rf_validation<-data.frame(

																																																															  				level=level,

																																																																							county=county,

																																																																							moYr=moYr,

																																																																											stationCount=length(anom_sp),

																																																																											krigeRFmmMax=krigeRFmax,

																																																																															krigeRFmmMin=krigeRFmin,

																																																																															negRFgridZero=negRFzeroGrid,

																																																																																			negRFpixCount=negRFpixCount,

																																																																																			gridPixCount=gridPixCount,

																																																																																							negRFloocv=as.logical(max(loocv_df$negRFzero)),

																																																																																							rsq_rf_mm=summary(lm(loocv_df$obs_rf~loocv_df$pred_rf))$r.squared,

																																																																																											rmse_rf_mm=rmse(loocv_df$obs_rf,loocv_df$pred_rf),

																																																																																											mae_rf_mm=mae(loocv_df$obs_rf,loocv_df$pred_rf),

																																																																																															bias_rf_mm=bias(loocv_df$obs_rf,loocv_df$pred_rf),

																																																																																															rsq_rf_anom=summary(lm(loocv_df$obs_anom~loocv_df$pred_anom))$r.squared,

																																																																																																			rmse_rf_anom=rmse(loocv_df$obs_anom,loocv_df$pred_anom),

																																																																																																			mae_rf_anom=mae(loocv_df$obs_anom,loocv_df$pred_anom),

																																																																																																							bias_rf_anom=bias(loocv_df$obs_anom,loocv_df$pred_anom)

																																																																																																							)



																																																												rf_fix_df<-cbind(names(varioList)[v],uglyClass,rf_validation)

																																																															names(rf_fix_df)[1]<-"vario"

																																																															rf_fix_df<-rf_fix_df[,c("county","moYr","vario","level","Class","probUgly","krigeRFmmMax","krigeRFmmMin","negRFgridZero","negRFpixCount","gridPixCount",

																																																																									"stationCount","negRFloocv","rsq_rf_mm","rmse_rf_mm","mae_rf_mm","bias_rf_mm","rsq_rf_anom","rmse_rf_anom","mae_rf_anom","bias_rf_anom")]#reorder cols

																																																																		#add to list to save

																																																																		blankList<-list()

																																																																					blankList[[paste0(county,"_anom_sp")]]<-anom_sp

																																																																					blankList[["variogram"]]<-vario

																																																																								blankList[["krigeObj"]]<-krigeObj

																																																																								blankList[["rf_mm_ras"]]<-rf_mm_ras

																																																																											blankList[["rf_mm_SE_ras"]]<-rf_mm_SE_ras

																																																																											blankList[["rf_anom_ras"]]<-rf_anom_ras

																																																																														blankList[["rf_anom_SE_ras"]]<-rf_anom_SE_ras

																																																																														blankList[["loocv_df"]]<-loocv_df

																																																																																	blankList[["map_validation_df"]]<-rf_fix_df

																																																																																	runList[[names(varioList)[v]]]<-blankList

																																																																																				#print(runList)

																																																																																				rf_fix_df_all<-rbind(rf_fix_df_all,rf_fix_df)

																																																																																				}#end vario loop

														runList[[paste0("rf_map_check_",level)]]<-rf_fix_df_all

														return(runList)

																}#end uglyMapFix function



#vario_name_func function

vario_name_func<-function(map_check,noFix=as.logical(0),freePriority=as.logical(0)){

	    rescaleRF <- function(x,rfmax) {100-(round(x,0)/round(rfmax,0) * 100)} #rescale rf mm metrics by max krige rf

    #noFix condition = FALSE default

    if(!noFix){

	         map_check_sub<-map_check[map_check$Class=="notUgly",] #subset all 'notUgly' maps

         if(nrow(map_check_sub)==0){return(as.character(map_check_sub$vario))}else{

		       if(as.logical(max(map_check_sub$vario=="vari_free_all")) & freePriority){ #if "free all" is not ugly choose free all map

			       	 best_vario_name<-as.character("vari_free_all")

	 	 }else{ #"free all" map isUgly rescale rmse & mae and combind with rsq to calc best rescaled median value map 

			 	 map_check_sub$rmseRe<-rescaleRF(map_check_sub$rmse_rf_mm,rfmax=map_check_sub$krigeRFmmMax)#rescale rmse rf based on highest krige rf

		 	 map_check_sub$maeRe<-rescaleRF(map_check_sub$mae_rf_mm,rfmax=map_check_sub$krigeRFmmMax)#rescale mae rf based on highest krige rf

			 	 map_check_sub$rsqRe<-round(map_check_sub$rsq_rf_mm,2)*100 #round and rescale rf 0-100

			 	 map_check_sub$probUglyRe<-round(1-map_check_sub$probUgly,2)*100 #reverse, round and rescale rf 0-100

				 	 map_check_sub$medianRe <- as.numeric(apply(map_check_sub[,c("rsqRe","maeRe","rmseRe","probUglyRe")], 1, median, na.rm = T))

				 	 map_check_sub<-map_check_sub[map_check_sub$medianRe==max(map_check_sub$medianRe),]#get highest median err/acc vario WITH lowest "probUgly"

					 	 best_vario_name<-as.character(map_check_sub$vario)

					 	 }

	    	 if(length(best_vario_name)>1){ #if more than 1 vario is tied for top median error/acc then use highest rsq vario within 0.01

			 		 map_check_sub<-map_check_sub[round(map_check_sub$rsq_rf_mm,2)==round(max(map_check_sub$rsq_rf_mm),2),]#get highest rsq rf mm vario within 0.01

		 	 	 best_vario_name<-as.character(map_check_sub$vario)

				 	 }

		    	 if(length(best_vario_name)>1){ #if STILL more than 1 vario is tied for best then use LOWEST mae within a mm

				 		 map_check_sub<-map_check_sub[round(map_check_sub$mae_rf_mm,0)==round(max(map_check_sub$mae_rf_mm),0),]#get lowest mae rf mm vario within 1 mm

			 	 	 best_vario_name<-as.character(map_check_sub$vario)

					 	 }

			    	 if(length(best_vario_name)>1){ #if STILL more than 1 vario is tied for best then use LOWEST rmse within a mm

					 		 map_check_sub<-map_check_sub[round(map_check_sub$rmse_rf_mm,0)==round(max(map_check_sub$rmse_rf_mm),0),]#get lowest rmse rf mm vario within 1 mm

				 		 best_vario_name<-as.character(map_check_sub$vario)

						 	 }

				 	 if(length(best_vario_name)>1){ #if STILL more than 1 varios for best subset a random vario from both remaing varios...

						 		 set.seed(22)#reproducable randomness

					 		 map_check_sub<-map_check_sub[sample(nrow(map_check_sub),1),] #if not then randomly choose from remaing varios

							 		 best_vario_name<-as.character(map_check_sub$vario)

							 		}

					 	 return(best_vario_name)

						 	 }}#end if default function setting

        

    	#noFix condition = TRUE

    	if(noFix){

			map_check_sub<-map_check

		 #all maps ARE "isUgly" rescale rmse & mae and combind with rsq and reversed probUgly to calc best case rescaled median value map 

		 map_check_sub$rmseRe<-rescaleRF(map_check_sub$rmse_rf_mm,rfmax=map_check_sub$krigeRFmmMax)#rescale rmse rf based on highest krige rf

		 	 map_check_sub$maeRe<-rescaleRF(map_check_sub$mae_rf_mm,rfmax=map_check_sub$krigeRFmmMax)#rescale mae rf based on highest krige rf

		 	 map_check_sub$rsqRe<-round(map_check_sub$rsq_rf_mm,2)*100 #round and rescale rf 0-100

			 	 map_check_sub$probUglyRe<-round(1-map_check_sub$probUgly,2)*100 #reverse, round and rescale rf 0-100

			 	 map_check_sub$medianRe <- as.numeric(apply(map_check_sub[,c("rsqRe","maeRe","rmseRe","probUglyRe")], 1, median, na.rm = T))

				 	 map_check_sub<-map_check_sub[map_check_sub$medianRe==max(map_check_sub$medianRe),]#get highest median err/acc vario WITH lowest "probUgly"

				 	 best_vario_name<-as.character(map_check_sub$vario)

					 	 #check if only 1 vario

					    	 if(length(best_vario_name)>1){ #if more than 1 vario is tied for top median error/acc then use lowest probUgly vario within 0.01

							 		 map_check_sub<-map_check_sub[round(map_check_sub$probUgly,2)==round(min(map_check_sub$probUgly),2),]#get lowest probUgly vario within 0.01

					 		 best_vario_name<-as.character(map_check_sub$vario)

							 	 }

					    	 if(length(best_vario_name)>1){ #if STILL more than 1 vario is tied for top median error/acc then use highest rsq vario within 0.01

							 		 map_check_sub<-map_check_sub[round(map_check_sub$rsq_rf_mm,2)==round(max(map_check_sub$rsq_rf_mm),2),]#get highest rsq rf mm vario within 0.01

						 		 best_vario_name<-as.character(map_check_sub$vario)

								 	 }

						    	 if(length(best_vario_name)>1){ #if STILL more than 1 vario is tied for best then use LOWEST mae within a mm

								 		 map_check_sub<-map_check_sub[round(map_check_sub$mae_rf_mm,0)==round(max(map_check_sub$mae_rf_mm),0),]#get lowest mae rf mm vario within 1 mm

							 		 best_vario_name<-as.character(map_check_sub$vario)

									 	 }

							    	 if(length(best_vario_name)>1){ #if STILL more than 1 vario is tied for best then use LOWEST rmse within a mm

									 		 map_check_sub<-map_check_sub[round(map_check_sub$rmse_rf_mm,0)==round(max(map_check_sub$rmse_rf_mm),0),]#get lowest rmse rf mm vario within 1 mm

								 		 best_vario_name<-as.character(map_check_sub$vario)

										 	 }

								 	 if(length(best_vario_name)>1){ #if STILL more than 1 varios for best subset a random vario from both remaing varios...

										 		set.seed(22)#reproducable randomness

									 		map_check_sub<-map_check_sub[sample(nrow(map_check_sub),1),] #if not then randomly choose from remaing varios

													best_vario_name<-as.character(map_check_sub$vario)

													}

									 	 return(best_vario_name)

										 	 }#end noFix setting

    	}#end all vario_name_func function



##multiplot all kriges





multiKrigePlot<-function(kriges,map_check_all,best_vario_name){



	 require(raster)

 getMaxVal<-function(ras){ifelse(is.null(ras),0,maxValue(ras))}



  var_level<-ifelse(max(unique(map_check_all$level)=="high"),"high",

		    		ifelse(max(unique(map_check_all$level)=="med"),"med",

					   			if(max(unique(map_check_all$level)=="low")){"low"}))



  ###IF bi...

  if(as.character(map_check_all$county[1])=="bi"){

	   

	   #define kriges and rasters

	   bi_kriges<-kriges

   bi_free_all_rf<-bi_kriges$vari_free_all$rf_mm_ras

    bi_fix_nug_rf<-bi_kriges[[paste0("vari_fix_nug_",var_level)]]$rf_mm_ras

    bi_fix_rng_rf<-bi_kriges[[paste0("vari_fix_rng_",var_level)]]$rf_mm_ras

     bi_free_sill_rf<-bi_kriges[[paste0("vari_free_sil_",var_level)]]$rf_mm_ras

     

     #re-make map check subset

     bi_map_check<-rbind(map_check_all[map_check_all$vario=="vari_free_all",],map_check_all[map_check_all$vario!="vari_free_all" & map_check_all$level==var_level,])

      bi_map_check$best<-bi_map_check$vario==best_vario_name

      

      #get max rf from all 4 rf rasters

      max_mm_bi<-max(c(getMaxVal(bi_free_all_rf),getMaxVal(bi_fix_nug_rf),getMaxVal(bi_fix_rng_rf),getMaxVal(bi_free_sill_rf)))

       

       #MAKE PLOTS 

        par(mfrow=c(2,2),mai = c(0.25,0.25,0.25,0.25))

        if(!is.null(bi_free_all_rf)){

		   plot(bi_free_all_rf,col=rainbow(300,end=0.8),zlim=c(0,max_mm_bi),axes=FALSE, box=FALSE,legend=FALSE,

			 	main="vari_free_all") 

	   text(x = extent(bi_free_all_rf)[2]-0.3,y = extent(bi_free_all_rf)[3]+0.2, labels = 

				paste0("Class: ",bi_map_check[bi_map_check$vario=="vari_free_all",]$Class,

				       			 "\n probUgly: ",round(bi_map_check[bi_map_check$vario=="vari_free_all",]$probUgly*100,1),

							 			 "%\n Rsq: ",round(bi_map_check[bi_map_check$vario=="vari_free_all",]$rsq_rf_mm,3),

							 			 "\n RMSE: ",round(bi_map_check[bi_map_check$vario=="vari_free_all",]$rmse_rf_mm,1),

										 			 "\n MAE: ",round(bi_map_check[bi_map_check$vario=="vari_free_all",]$mae_rf_mm,1)))

	      if(bi_map_check[bi_map_check$vario=="vari_free_all","best"]){

		      	text(x = extent(bi_free_all_rf)[2]-0.3,y = extent(bi_free_all_rf)[4]-0.13, labels = "*BEST*", col="blue")

	       	}

	      if(bi_map_check[bi_map_check$vario=="vari_free_all",]$negRFgrid){

		      	text(x = extent(bi_free_all_rf)[2]-0.3,y = extent(bi_free_all_rf)[4]-0.22, labels = "*Neg RF Fix*", col="red")

	       	}

	        }else{

			   plot(0,type='n',axes=FALSE,ann=FALSE,xlim=c(0,0),ylim=c(0,0))

		   text(x = 0,y = 0, labels = "missing", col="red")

		     }



	  if(!is.null(bi_fix_nug_rf)){

		     plot(bi_fix_nug_rf,col=rainbow(300,end=0.8),zlim=c(0,max_mm_bi),axes=FALSE, box=FALSE,legend=FALSE,

			   	main=paste0("vari_fix_nug_",var_level)) 

	     text(x = extent(bi_fix_nug_rf)[2]-0.3,y = extent(bi_fix_nug_rf)[3]+0.2, labels = 

		  		paste0("Class: ",bi_map_check[bi_map_check$vario==paste0("vari_fix_nug_",var_level),]$Class,

					 			 "\n probUgly: ",round(bi_map_check[bi_map_check$vario==paste0("vari_fix_nug_",var_level),]$probUgly*100,1),

								 			 "%\n Rsq: ",round(bi_map_check[bi_map_check$vario==paste0("vari_fix_nug_",var_level),]$rsq_rf_mm,3),

								 			 "\n RMSE: ",round(bi_map_check[bi_map_check$vario==paste0("vari_fix_nug_",var_level),]$rmse_rf_mm,1),

											 			 "\n MAE: ",round(bi_map_check[bi_map_check$vario==paste0("vari_fix_nug_",var_level),]$mae_rf_mm,1)))

	        if(bi_map_check[bi_map_check$vario==paste0("vari_fix_nug_",var_level),"best"]){

				text(x = extent(bi_free_all_rf)[2]-0.3,y = extent(bi_free_all_rf)[4]-0.13, labels = "*BEST*", col="blue")

		 	}

	        if(bi_map_check[bi_map_check$vario==paste0("vari_fix_nug_",var_level),]$negRFgrid){

				text(x = extent(bi_free_all_rf)[2]-0.3,y = extent(bi_free_all_rf)[4]-0.22, labels = "*Neg RF Fix*", col="red")

		 	}

		  }else{

			     plot(0,type='n',axes=FALSE,ann=FALSE,xlim=c(0,0),ylim=c(0,0))

		     text(x = 0,y = 0, labels = "missing", col="red")

		       }



	  if(!is.null(bi_fix_rng_rf)){

		    plot(bi_fix_rng_rf,col=rainbow(300,end=0.8),zlim=c(0,max_mm_bi),axes=FALSE, box=FALSE,legend=FALSE,

			  	main=paste0("vari_fix_rng_",var_level))

	    text(x = extent(bi_fix_rng_rf)[2]-0.3,y = extent(bi_fix_rng_rf)[3]+0.2, labels = 

		 		paste0("Class: ",bi_map_check[bi_map_check$vario==paste0("vari_fix_rng_",var_level),]$Class,

								 "\n probUgly: ",round(bi_map_check[bi_map_check$vario==paste0("vari_fix_rng_",var_level),]$probUgly*100,1),

								 			 "%\n Rsq: ",round(bi_map_check[bi_map_check$vario==paste0("vari_fix_rng_",var_level),]$rsq_rf_mm,3),

								 			 "\n RMSE: ",round(bi_map_check[bi_map_check$vario==paste0("vari_fix_rng_",var_level),]$rmse_rf_mm,1),

											 			 "\n mae: ",round(bi_map_check[bi_map_check$vario==paste0("vari_fix_rng_",var_level),]$mae_rf_mm,1)))

	      if(bi_map_check[bi_map_check$vario==paste0("vari_fix_rng_",var_level),"best"]){

		      	text(x = extent(bi_free_all_rf)[2]-0.3,y = extent(bi_free_all_rf)[4]-0.13, labels = "*BEST*", col="blue")

	       	}

	       if(bi_map_check[bi_map_check$vario==paste0("vari_fix_rng_",var_level),]$negRFgrid){

		       	text(x = extent(bi_free_all_rf)[2]-0.3,y = extent(bi_free_all_rf)[4]-0.22, labels = "*Neg RF Fix*", col="red")

	        	}

	        }else{

			   plot(0,type='n',axes=FALSE,ann=FALSE,xlim=c(0,0),ylim=c(0,0))

		   text(x = 0,y = 0, labels = "missing", col="red")

		     }



	    if(!is.null(bi_free_sill_rf)){

		       plot(bi_free_sill_rf,col=rainbow(300,end=0.8),zlim=c(0,max_mm_bi),axes=FALSE, box=FALSE,legend=FALSE,

			     	main=paste0("vari_free_sil_",var_level))

	       text(x = extent(bi_free_sill_rf)[2]-0.3,y = extent(bi_free_sill_rf)[3]+0.2, labels =

		    			paste0("Class: ",bi_map_check[bi_map_check$vario==paste0("vari_free_sil_",var_level),]$Class,

						   				 "\n probUgly: ",round(bi_map_check[bi_map_check$vario==paste0("vari_free_sil_",var_level),]$probUgly*100,1),

										 				 "%\n Rsq: ",round(bi_map_check[bi_map_check$vario==paste0("vari_free_sil_",var_level),]$rsq_rf_mm,3),

										 				 "\n RMSE: ",round(bi_map_check[bi_map_check$vario==paste0("vari_free_sil_",var_level),]$rmse_rf_mm,1),

														 				 "\n MAE: ",round(bi_map_check[bi_map_check$vario==paste0("vari_free_sil_",var_level),]$mae_rf_mm,1)))

	          if(bi_map_check[bi_map_check$vario==paste0("vari_free_sil_",var_level),"best"]){

			  	text(x = extent(bi_free_all_rf)[2]-0.3,y = extent(bi_free_all_rf)[4]-0.13, labels = "*BEST*", col="blue")

		   	}

	          if(bi_map_check[bi_map_check$vario==paste0("vari_free_sil_",var_level),]$negRFgrid){

			  	text(x = extent(bi_free_all_rf)[2]-0.3,y = extent(bi_free_all_rf)[4]-0.22, labels = "*Neg RF Fix*", col="red")

		   	}

		    }else{

			       plot(0,type='n',axes=FALSE,ann=FALSE,xlim=c(0,0),ylim=c(0,0))

		       text(x = 0,y = 0, labels = "missing", col="red")

		         }

	   }#if bi end plot





   ###IF mn...

   if(as.character(map_check_all$county[1])=="mn"){



	    #define kriges and rasters

	    mn_kriges<-kriges

    mn_free_all_rf<-mn_kriges$vari_free_all$rf_mm_ras

     mn_fix_nug_rf<-mn_kriges[[paste0("vari_fix_nug_",var_level)]]$rf_mm_ras

     mn_fix_rng_rf<-mn_kriges[[paste0("vari_fix_rng_",var_level)]]$rf_mm_ras

      mn_free_sill_rf<-mn_kriges[[paste0("vari_free_sil_",var_level)]]$rf_mm_ras

      

      #re-make map check subset

      mn_map_check<-rbind(map_check_all[map_check_all$vario=="vari_free_all",],map_check_all[map_check_all$vario!="vari_free_all" & map_check_all$level==var_level,])

       mn_map_check$best<-mn_map_check$vario==best_vario_name

       

       #get max rf from all 4 rf rasters

       max_mm_mn<-max(c(getMaxVal(mn_free_all_rf),getMaxVal(mn_fix_nug_rf),getMaxVal(mn_fix_rng_rf),getMaxVal(mn_free_sill_rf)))

        

        #MAKE PLOTS

         par(mfrow=c(2,2),mai = c(0.25,0.25,0.25,0.25))

         if(!is.null(mn_free_all_rf)){

		    plot(mn_free_all_rf,col=rainbow(300,end=0.8),zlim=c(0,max_mm_mn),axes=FALSE, box=FALSE,legend=FALSE,

			  	main="vari_free_all") 

	    text(x = extent(mn_free_all_rf)[1]+0.3,y = extent(mn_free_all_rf)[3]+0.2, labels = 

		 		paste0("Class: ",mn_map_check[mn_map_check$vario=="vari_free_all",]$Class,

								 "\n probUgly: ",round(mn_map_check[mn_map_check$vario=="vari_free_all",]$probUgly*100,1),

								 			 "%\n Rsq: ",round(mn_map_check[mn_map_check$vario=="vari_free_all",]$rsq_rf_mm,3),

								 			 "\n RMSE: ",round(mn_map_check[mn_map_check$vario=="vari_free_all",]$rmse_rf_mm,1),

											 			 "\n MAE: ",round(mn_map_check[mn_map_check$vario=="vari_free_all",]$mae_rf_mm,1)))



	       if(mn_map_check[mn_map_check$vario=="vari_free_all","best"]){

		       	text(x = extent(mn_free_all_rf)[2]-0.3,y = extent(mn_free_all_rf)[4]-0.13, labels = "*BEST*", col="blue")

	        	}

	       if(mn_map_check[mn_map_check$vario=="vari_free_all",]$negRFgrid){

		       	text(x = extent(mn_free_all_rf)[2]-0.3,y = extent(mn_free_all_rf)[4]-0.22, labels = "*NEG RF FIX*", col="red")

	        	}

	         }else{

			    plot(0,type='n',axes=FALSE,ann=FALSE,xlim=c(0,0),ylim=c(0,0))

		    text(x = 0,y = 0, labels = "missing", col="red")

		      }



	   if(!is.null(mn_fix_nug_rf)){

		      plot(mn_fix_nug_rf,col=rainbow(300,end=0.8),zlim=c(0,max_mm_mn),axes=FALSE, box=FALSE,legend=FALSE,

			    	main=paste0("vari_fix_nug_",var_level)) 

	       text(x = extent(mn_free_all_rf)[1]+0.3,y = extent(mn_free_all_rf)[3]+0.2, labels = 

		    		paste0("Class: ",mn_map_check[mn_map_check$vario==paste0("vari_fix_nug_",var_level),]$Class,

					   			 "\n probUgly: ",round(mn_map_check[mn_map_check$vario==paste0("vari_fix_nug_",var_level),]$probUgly*100,1),

								 			 "%\n Rsq: ",round(mn_map_check[mn_map_check$vario==paste0("vari_fix_nug_",var_level),]$rsq_rf_mm,3),

								 			 "\n RMSE: ",round(mn_map_check[mn_map_check$vario==paste0("vari_fix_nug_",var_level),]$rmse_rf_mm,1),

											 			 "\n MAE: ",round(mn_map_check[mn_map_check$vario==paste0("vari_fix_nug_",var_level),]$mae_rf_mm,1)))

	           if(mn_map_check[mn_map_check$vario==paste0("vari_fix_nug_",var_level),"best"]){

			   	text(x = extent(mn_free_all_rf)[2]-0.3,y = extent(mn_free_all_rf)[4]-0.13, labels = "*BEST*", col="blue")

		    	}

	           if(mn_map_check[mn_map_check$vario==paste0("vari_fix_nug_",var_level),]$negRFgrid){

			    	text(x = extent(mn_fix_nug_rf)[2]-0.3,y = extent(mn_fix_nug_rf)[4]-0.22, labels = "*NEG RF FIX*", col="red")

		    	}

		      }else{

			          plot(0,type='n',axes=FALSE,ann=FALSE,xlim=c(0,0),ylim=c(0,0))

		          text(x = 0,y = 0, labels = "missing", col="red")

			     }



	   if(!is.null(mn_fix_rng_rf)){

		      plot(mn_fix_rng_rf,col=rainbow(300,end=0.8),zlim=c(0,max_mm_mn),axes=FALSE, box=FALSE,legend=FALSE,

			    	main=paste0("vari_fix_rng_",var_level))

	      text(x = extent(mn_free_all_rf)[1]+0.3,y = extent(mn_free_all_rf)[3]+0.2, labels = 

		   		paste0("Class: ",mn_map_check[mn_map_check$vario==paste0("vari_fix_rng_",var_level),]$Class,

					  			 "\n probUgly: ",round(mn_map_check[mn_map_check$vario==paste0("vari_fix_rng_",var_level),]$probUgly*100,1),

								 			 "%\n Rsq: ",round(mn_map_check[mn_map_check$vario==paste0("vari_fix_rng_",var_level),]$rsq_rf_mm,3),

								 			 "\n RMSE: ",round(mn_map_check[mn_map_check$vario==paste0("vari_fix_rng_",var_level),]$rmse_rf_mm,1),

											 			 "\n MAE: ",round(mn_map_check[mn_map_check$vario==paste0("vari_fix_rng_",var_level),]$mae_rf_mm,1)))

	         if(mn_map_check[mn_map_check$vario==paste0("vari_fix_rng_",var_level),"best"]){

			 	text(x = extent(mn_free_all_rf)[2]-0.3,y = extent(mn_free_all_rf)[4]-0.13, labels = "*BEST*", col="blue")

		  	}

	         if(mn_map_check[mn_map_check$vario==paste0("vari_fix_rng_",var_level),]$negRFgrid){

			 	text(x = extent(mn_fix_rng_rf)[2]-0.3,y = extent(mn_fix_rng_rf)[4]-0.22, labels = "*NEG RF FIX*", col="red")

		    }

		   }else{

			      plot(0,type='n',axes=FALSE,ann=FALSE,xlim=c(0,0),ylim=c(0,0))

		      text(x = 0,y = 0, labels = "missing", col="red")

		        }

	     

	    if(!is.null(mn_free_sill_rf)){

		      plot(mn_free_sill_rf,col=rainbow(300,end=0.8),zlim=c(0,max_mm_mn),axes=FALSE, box=FALSE,legend=FALSE,

			    	main=paste0("vari_free_sil_",var_level))

	      text(x = extent(mn_free_all_rf)[1]+0.3,y = extent(mn_free_all_rf)[3]+0.2,, labels =

		   			paste0("Class: ",mn_map_check[mn_map_check$vario==paste0("vari_free_sil_",var_level),]$Class,

						  				 "\n probUgly: ",round(mn_map_check[mn_map_check$vario==paste0("vari_free_sil_",var_level),]$probUgly*100,1),

										 				 "%\n Rsq: ",round(mn_map_check[mn_map_check$vario==paste0("vari_free_sil_",var_level),]$rsq_rf_mm,3),

										 				 "\n RMSE: ",round(mn_map_check[mn_map_check$vario==paste0("vari_free_sil_",var_level),]$rmse_rf_mm,1),

														 			 	 "\n MAE: ",round(mn_map_check[mn_map_check$vario==paste0("vari_free_sil_",var_level),]$mae_rf_mm,1)))

	        if(mn_map_check[mn_map_check$vario==paste0("vari_free_sil_",var_level),"best"]){

				text(x = extent(mn_free_all_rf)[2]-0.3,y = extent(mn_free_all_rf)[4]-0.13, labels = "*BEST*", col="blue")

		 	}

	        if(mn_map_check[mn_map_check$vario==paste0("vari_free_sil_",var_level),]$negRFgrid){

				text(x = extent(mn_free_sill_rf)[2]-0.3,y = extent(mn_free_sill_rf)[4]-0.22, labels = "*NEG RF FIX*", col="red")

		    }

		   }else{

			       plot(0,type='n',axes=FALSE,ann=FALSE,xlim=c(0,0),ylim=c(0,0))

		       text(x = 0,y = 0, labels = "missing", col="red")

		          }

	     }#if mn end plot



  ###IF oa...

   if(as.character(map_check_all$county[1])=="oa"){

	    

	    #define kriges and rasters

	    oa_kriges<-kriges

    oa_free_all_rf<-oa_kriges$vari_free_all$rf_mm_ras

     oa_fix_nug_rf<-oa_kriges[[paste0("vari_fix_nug_",var_level)]]$rf_mm_ras

     oa_fix_rng_rf<-oa_kriges[[paste0("vari_fix_rng_",var_level)]]$rf_mm_ras

      oa_free_sill_rf<-oa_kriges[[paste0("vari_free_sil_",var_level)]]$rf_mm_ras

      

      #re-make map check subset

      oa_map_check<-rbind(map_check_all[map_check_all$vario=="vari_free_all",],map_check_all[map_check_all$vario!="vari_free_all" & map_check_all$level==var_level,])

       oa_map_check$best<-oa_map_check$vario==best_vario_name

       

       #get max rf from all 4 rf rasters

       max_mm_oa<-max(c(getMaxVal(oa_free_all_rf),getMaxVal(oa_fix_nug_rf),getMaxVal(oa_fix_rng_rf),getMaxVal(oa_free_sill_rf)))

        

        #MAKE PLOTS

         par(mfrow=c(2,2),mai = c(0.25,0.25,0.25,0.25))

        if(!is.null(oa_free_all_rf)){

		   plot(oa_free_all_rf,col=rainbow(300,end=0.8),zlim=c(0,max_mm_oa),axes=FALSE, box=FALSE,legend=FALSE,

			 	main="vari_free_all") 

	   text(x = extent(oa_free_all_rf)[2]-0.13,y = extent(oa_free_all_rf)[4]- 0.09, labels = 

		 		paste0("Class: ",oa_map_check[oa_map_check$vario=="vari_free_all",]$Class,

					 			 "\n probUgly: ",round(oa_map_check[oa_map_check$vario=="vari_free_all",]$probUgly*100,1),

								 			 "%\n Rsq: ",round(oa_map_check[oa_map_check$vario=="vari_free_all",]$rsq_rf_mm,3),

								 			 "\n RMSE: ",round(oa_map_check[oa_map_check$vario=="vari_free_all",]$rmse_rf_mm,1),

											 			 "\n MAE: ",round(oa_map_check[oa_map_check$vario=="vari_free_all",]$mae_rf_mm,1)))

	      if(oa_map_check[oa_map_check$vario=="vari_free_all","best"]){

		      	text(x = extent(oa_free_all_rf)[1]+0.15,y = extent(oa_free_all_rf)[4]-0.07, labels = "*BEST*", col="blue")

	       	}

	      if(oa_map_check[oa_map_check$vario=="vari_free_all",]$negRFgrid){

		      	text(x = extent(oa_free_all_rf)[1]+0.15,y = extent(oa_free_all_rf)[4]-0.11, labels = "*NEG RF FIX*", col="red")

	         }

	         }else{

			     plot(0,type='n',axes=FALSE,ann=FALSE,xlim=c(0,0),ylim=c(0,0))

		     text(x = 0,y = 0, labels = "missing", col="red")

		        }



	 if(!is.null(oa_fix_nug_rf)){

		   plot(oa_fix_nug_rf,col=rainbow(300,end=0.8),zlim=c(0,max_mm_oa),axes=FALSE, box=FALSE,legend=FALSE,

			 	main=paste0("vari_fix_nug_",var_level)) 

	   text(x = extent(oa_free_all_rf)[2]-0.13,y = extent(oa_free_all_rf)[4]- 0.09, labels = 

				paste0("Class: ",oa_map_check[oa_map_check$vario==paste0("vari_fix_nug_",var_level),]$Class,

				       			 "\n probUgly: ",round(oa_map_check[oa_map_check$vario==paste0("vari_fix_nug_",var_level),]$probUgly*100,1),

							 			 "%\n Rsq: ",round(oa_map_check[oa_map_check$vario==paste0("vari_fix_nug_",var_level),]$rsq_rf_mm,3),

							 			 "\n RMSE: ",round(oa_map_check[oa_map_check$vario==paste0("vari_fix_nug_",var_level),]$rmse_rf_mm,1),

										 			 "\n MAE: ",round(oa_map_check[oa_map_check$vario==paste0("vari_fix_nug_",var_level),]$mae_rf_mm,1)))



	     if(oa_map_check[oa_map_check$vario==paste0("vari_fix_nug_",var_level),"best"]){

		     	text(x = extent(oa_free_all_rf)[1]+0.15,y = extent(oa_free_all_rf)[4]-0.07, labels = "*BEST*", col="blue")

	      	}

	     if(oa_map_check[oa_map_check$vario==paste0("vari_fix_nug_",var_level),]$negRFgrid){

		     	text(x = extent(oa_free_all_rf)[1]+0.15,y = extent(oa_free_all_rf)[4]-0.11, labels = "*NEG RF FIX*", col="red")

	      	}

	        }else{

			    plot(0,type='n',axes=FALSE,ann=FALSE,xlim=c(0,0),ylim=c(0,0))

		    text(x = 0,y = 0, labels = "missing", col="red")

		       }



	 if(!is.null(oa_fix_rng_rf)){

		   plot(oa_fix_rng_rf,col=rainbow(300,end=0.8),zlim=c(0,max_mm_oa),axes=FALSE, box=FALSE,legend=FALSE,

			 	main=paste0("vari_fix_rng_",var_level))

	   text(x = extent(oa_free_all_rf)[2]-0.13,y = extent(oa_free_all_rf)[4]- 0.09, labels = 

				paste0("Class: ",oa_map_check[oa_map_check$vario==paste0("vari_fix_rng_",var_level),]$Class,

				       			 "\n probUgly: ",round(oa_map_check[oa_map_check$vario==paste0("vari_fix_rng_",var_level),]$probUgly*100,1),

							 			 "%\n Rsq: ",round(oa_map_check[oa_map_check$vario==paste0("vari_fix_rng_",var_level),]$rsq_rf_mm,3),

							 			 "\n RMSE: ",round(oa_map_check[oa_map_check$vario==paste0("vari_fix_rng_",var_level),]$rmse_rf_mm,1),

										 			 "\n MAE: ",round(oa_map_check[oa_map_check$vario==paste0("vari_fix_rng_",var_level),]$mae_rf_mm,1)))

	     if(oa_map_check[oa_map_check$vario==paste0("vari_fix_rng_",var_level),"best"]){

		     	text(x = extent(oa_free_all_rf)[1]+0.15,y = extent(oa_free_all_rf)[4]-0.07, labels = "*BEST*", col="blue")

	      	}

	     if(oa_map_check[oa_map_check$vario==paste0("vari_fix_rng_",var_level),]$negRFgrid){

		     	text(x = extent(oa_free_all_rf)[1]+0.15,y = extent(oa_free_all_rf)[4]-0.11, labels = "*NEG RF FIX*", col="red")

	         }

	         }else{

			     plot(0,type='n',axes=FALSE,ann=FALSE,xlim=c(0,0),ylim=c(0,0))

		     text(x = 0,y = 0, labels = "missing", col="red")

		        }



	  if(!is.null(oa_free_sill_rf)){

		    plot(oa_free_sill_rf,col=rainbow(300,end=0.8),zlim=c(0,max_mm_oa),axes=FALSE, box=FALSE,legend=FALSE,

			  	main=paste0("vari_free_sil_",var_level))

	    text(x = extent(oa_free_all_rf)[2]-0.13,y = extent(oa_free_all_rf)[4]- 0.09, labels =

		 			paste0("Class: ",oa_map_check[oa_map_check$vario==paste0("vari_free_sil_",var_level),]$Class,

										 "\n probUgly: ",round(oa_map_check[oa_map_check$vario==paste0("vari_free_sil_",var_level),]$probUgly*100,1),

										 				 "%\n Rsq: ",round(oa_map_check[oa_map_check$vario==paste0("vari_free_sil_",var_level),]$rsq_rf_mm,3),

										 				 "\n RMSE: ",round(oa_map_check[oa_map_check$vario==paste0("vari_free_sil_",var_level),]$rmse_rf_mm,1),

														 			 	 "\n MAE: ",round(oa_map_check[oa_map_check$vario==paste0("vari_free_sil_",var_level),]$mae_rf_mm,1)))

	      if(oa_map_check[oa_map_check$vario==paste0("vari_free_sil_",var_level),"best"]){

		      	text(x = extent(oa_free_all_rf)[1]+0.15,y = extent(oa_free_all_rf)[4]-0.07, labels = "*BEST*", col="blue")

	       	}

	      if(oa_map_check[oa_map_check$vario==paste0("vari_free_sil_",var_level),]$negRFgrid){

		      	text(x = extent(oa_free_all_rf)[1]+0.15,y = extent(oa_free_all_rf)[4]-0.11, labels = "*NEG RF FIX*", col="red")

	         }

	         }else{

			     plot(0,type='n',axes=FALSE,ann=FALSE,xlim=c(0,0),ylim=c(0,0))

		     text(x = 0,y = 0, labels = "missing", col="red")

		        }

	   }#if oa end plot





    ###IF ka...

    if(as.character(map_check_all$county[1])=="ka"){



	     #define kriges and rasters

	     ka_kriges<-kriges

     ka_free_all_rf<-ka_kriges$vari_free_all$rf_mm_ras

      ka_fix_nug_rf<-ka_kriges[[paste0("vari_fix_nug_",var_level)]]$rf_mm_ras

      ka_fix_rng_rf<-ka_kriges[[paste0("vari_fix_rng_",var_level)]]$rf_mm_ras

       ka_free_sill_rf<-ka_kriges[[paste0("vari_free_sil_",var_level)]]$rf_mm_ras



       #re-make map check subset

       ka_map_check<-rbind(map_check_all[map_check_all$vario=="vari_free_all",],map_check_all[map_check_all$vario!="vari_free_all" & map_check_all$level==var_level,])

        ka_map_check$best<-ka_map_check$vario==best_vario_name

        

        #get max rf from all 4 rf rasters

        max_mm_ka<-max(c(getMaxVal(ka_free_all_rf),getMaxVal(ka_fix_nug_rf),getMaxVal(ka_fix_rng_rf),getMaxVal(ka_free_sill_rf)))

	 

	 #MAKE PLOTS

	  par(mfrow=c(2,2),mai = c(0.25,0.25,0.25,0.25))

	 if(!is.null(ka_free_all_rf)){

		    plot(ka_free_all_rf,col=rainbow(300,end=0.8),zlim=c(0,max_mm_ka),axes=FALSE, box=FALSE,legend=FALSE,

			 	main="vari_free_all") 

	    text(x = extent(ka_free_all_rf)[1]+0.08,y = extent(ka_free_all_rf)[3]+0.09, labels = 

		 		paste0("Class:",ka_map_check[ka_map_check$vario=="vari_free_all",]$Class,

								 "\n probUgly:",round(ka_map_check[ka_map_check$vario=="vari_free_all",]$probUgly*100,1),"%"))

	       text(x = extent(ka_free_all_rf)[1]+0.14,y = extent(ka_free_all_rf)[3]+0.045, labels =

		    			paste0("\n Rsq:",round(ka_map_check[ka_map_check$vario=="vari_free_all",]$rsq_rf_mm,3),

						   			 "\n RMSE:",round(ka_map_check[ka_map_check$vario=="vari_free_all",]$rmse_rf_mm,1),

									 			 "; MAE: ",round(ka_map_check[ka_map_check$vario=="vari_free_all",]$mae_rf_mm,1)))

	       if(ka_map_check[ka_map_check$vario=="vari_free_all","best"]){

		       	text(x = extent(ka_free_all_rf)[1]+0.11,y = extent(ka_free_all_rf)[4]-0.03, labels = "*BEST*", col="blue")

	        	}

	          if(ka_map_check[ka_map_check$vario=="vari_free_all",]$negRFgrid){

			  	text(x = extent(ka_free_all_rf)[1]+0.11,y = extent(ka_free_all_rf)[4]-0.07, labels = "*NEG RF \n FIX*", col="red")

		   	}

	          }else{

			      plot(0,type='n',axes=FALSE,ann=FALSE,xlim=c(0,0),ylim=c(0,0))

		      text(x = 0,y = 0, labels = "missing", col="red")

		         }



	  if(!is.null(ka_fix_nug_rf)){

		     plot(ka_fix_nug_rf,col=rainbow(300,end=0.8),zlim=c(0,max_mm_ka),axes=FALSE, box=FALSE,legend=FALSE,

			   	main=paste0("vari_fix_nug_",var_level)) 

	      text(x = extent(ka_free_all_rf)[1]+0.08,y = extent(ka_free_all_rf)[3]+0.09, labels = 

		   		paste0("Class:",ka_map_check[ka_map_check$vario==paste0("vari_fix_nug_",var_level),]$Class,

					  			 "\n probUgly:",round(ka_map_check[ka_map_check$vario==paste0("vari_fix_nug_",var_level),]$probUgly*100,1),"%"))

	          text(x = extent(ka_free_all_rf)[1]+0.14,y = extent(ka_free_all_rf)[3]+0.045, labels =

		       			paste0("\n Rsq:",round(ka_map_check[ka_map_check$vario==paste0("vari_fix_nug_",var_level),]$rsq_rf_mm,3),

						      			 "\n RMSE:",round(ka_map_check[ka_map_check$vario==paste0("vari_fix_nug_",var_level),]$rmse_rf_mm,1),

									 			 "; MAE: ",round(ka_map_check[ka_map_check$vario==paste0("vari_fix_nug_",var_level),]$mae_rf_mm,1)))

	          if(ka_map_check[ka_map_check$vario==paste0("vari_fix_nug_",var_level),"best"]){

			  	text(x = extent(ka_free_all_rf)[1]+0.11,y = extent(ka_free_all_rf)[4]-0.03, labels = "*BEST*", col="blue")

		   	}

		      if(ka_map_check[ka_map_check$vario==paste0("vari_fix_nug_",var_level),]$negRFgrid){

			      	text(x = extent(ka_free_all_rf)[1]+0.11,y = extent(ka_free_all_rf)[4]-0.07, labels = "*NEG RF \n FIX*", col="red")

		       	}

		     }else{

			         plot(0,type='n',axes=FALSE,ann=FALSE,xlim=c(0,0),ylim=c(0,0))

		         text(x = 0,y = 0, labels = "missing", col="red")

			    }



	  if(!is.null(ka_fix_rng_rf)){

		     plot(ka_fix_rng_rf,col=rainbow(300,end=0.8),zlim=c(0,max_mm_ka),axes=FALSE, box=FALSE,legend=FALSE,

			   	main=paste0("vari_fix_rng_",var_level))

	     text(x = extent(ka_free_all_rf)[1]+0.08,y = extent(ka_free_all_rf)[3]+0.09, labels = 

		  		paste0("Class:",ka_map_check[ka_map_check$vario==paste0("vari_fix_rng_",var_level),]$Class,

					 			 "\n probUgly:",round(ka_map_check[ka_map_check$vario==paste0("vari_fix_rng_",var_level),]$probUgly*100,1),"%"))

	        text(x = extent(ka_free_all_rf)[1]+0.14,y = extent(ka_free_all_rf)[3]+0.045, labels =

		     			paste0("\n Rsq:",round(ka_map_check[ka_map_check$vario==paste0("vari_fix_rng_",var_level),]$rsq_rf_mm,3),

						    			 "\n RMSE:",round(ka_map_check[ka_map_check$vario==paste0("vari_fix_rng_",var_level),]$rmse_rf_mm,1),

									 			 "; MAE: ",round(ka_map_check[ka_map_check$vario==paste0("vari_fix_rng_",var_level),]$mae_rf_mm,1)))

	         if(ka_map_check[ka_map_check$vario==paste0("vari_fix_rng_",var_level),"best"]){

			 	text(x = extent(ka_free_all_rf)[1]+0.11,y = extent(ka_free_all_rf)[4]-0.03, labels = "*BEST*", col="blue")

		  	}

		    if(ka_map_check[ka_map_check$vario==paste0("vari_fix_rng_",var_level),]$negRFgrid){

			    	text(x = extent(ka_free_all_rf)[1]+0.11,y = extent(ka_free_all_rf)[4]-0.07, labels = "*NEG RF \n FIX*", col="red")

		    	}

		    }else{

			        plot(0,type='n',axes=FALSE,ann=FALSE,xlim=c(0,0),ylim=c(0,0))

		        text(x = 0,y = 0, labels = "missing", col="red")

			   }



	   if(!is.null(ka_free_sill_rf)){

		      plot(ka_free_sill_rf,col=rainbow(300,end=0.8),zlim=c(0,max_mm_ka),axes=FALSE, box=FALSE,legend=FALSE,

			    	main=paste0("vari_free_sil_",var_level))

	      text(x = extent(ka_free_all_rf)[1]+0.08,y = extent(ka_free_all_rf)[3]+0.09, labels = 

		   		paste0("Class:",ka_map_check[ka_map_check$vario==paste0("vari_free_sil_",var_level),]$Class,

					  			 "\n probUgly: ",round(ka_map_check[ka_map_check$vario==paste0("vari_free_sil_",var_level),]$probUgly*100,1),"%"))

	         text(x = extent(ka_free_all_rf)[1]+0.14,y = extent(ka_free_all_rf)[3]+0.045, labels =

		      			paste0("\n Rsq:",round(ka_map_check[ka_map_check$vario==paste0("vari_free_sil_",var_level),]$rsq_rf_mm,3),

						     			 "\n RMSE:",round(ka_map_check[ka_map_check$vario==paste0("vari_free_sil_",var_level),]$rmse_rf_mm,1),

									 			 "; MAE: ",round(ka_map_check[ka_map_check$vario==paste0("vari_free_sil_",var_level),]$mae_rf_mm,1)))

	         if(ka_map_check[ka_map_check$vario==paste0("vari_free_sil_",var_level),"best"]){

			 	text(x = extent(ka_free_all_rf)[1]+0.11,y = extent(ka_free_all_rf)[4]-0.03, labels = "*BEST*", col="blue")

		  	}

		    if(ka_map_check[ka_map_check$vario==paste0("vari_free_sil_",var_level),]$negRFgrid){

			    	text(x = extent(ka_free_all_rf)[1]+0.11,y = extent(ka_free_all_rf)[4]-0.07, labels = "*NEG RF \n FIX*", col="red")

		    	}

		    }else{

			        plot(0,type='n',axes=FALSE,ann=FALSE,xlim=c(0,0),ylim=c(0,0))

		        text(x = 0,y = 0, labels = "missing", col="red")

			   }

	    }#if ka end plot

}#end multiKrigePlot function



###metamaker function 

metamaker<-function(map_validation_df,grid,filenames,datatype,mon_yr=data_mon_yr,statewide=F,loocv_df=NA){

	  #packages

		require(raster)

  #get dates

	dataStartDate<-as.Date(paste(substr(data_mon_yr,1,4),as.integer(substrRight(data_mon_yr,2)),"1",sep="-"))

		dataEndDate<-seq(as.Date(paste(substr(data_mon_yr,1,4),as.integer(substrRight(data_mon_yr,2)),"1",sep="-")),length=2,by="months")[2]-1

		dataconvert<-function(x){return(ifelse(is.numeric(x),as.character(round(x,5)),as.character(x)))}

			

			if(statewide==F){

					valid_meta<-data.frame(attribute=as.character(names(map_validation_df)),value=as.character(lapply(map_validation_df[1,], dataconvert)))



			co_statement<-as.character(paste("This",format(dataStartDate,"%B %Y"),"rainfall map of",

							 		ifelse(map_validation_df$county=="ka","Kauai County",

												 ifelse(map_validation_df$county=="oa","Honolulu County (Oahu)",

															 ifelse(map_validation_df$county=="mn","Maui County (Maui, Lanai, Molokai & Kahoolawe)",

																		 if(map_validation_df$county=="bi"){"Hawaii county"}))), 

							 		"is a high spatial resolution (~250m) gridded prediction of cumulative rainfall in millimeters from",

											format(dataStartDate,"%d %b %Y"),"to",format(dataEndDate,"%d %b %Y."),

											"This was produced using a climate-aided modified automatic kriging interpolation of a log transformed sum of a monthly rainfall anomaly ratio (observed mm / mean monthly mm) and a constant of 1 (k). This kriging process used", 

													map_validation_df$stationCount, "unique station locations within",

													ifelse(map_validation_df$county=="ka","Kauai County",

													       		 ifelse(map_validation_df$county=="oa","Honolulu County (Oahu)",

																       	 	ifelse(map_validation_df$county=="mn","Maui County (Maui, Lanai, Molokai & Kahoolawe)",

																			       	 	if(map_validation_df$county=="bi"){"Hawaii county"}))),

							 		"and their",format(dataStartDate,"%B %Y"), 

											"recorded and/or estimated rainfall (mm) totals. A leave one out cross validation (LOOCV) of the station data used in this map produced an R-squared of:",

											round(map_validation_df$rsq_rf_mm,3),",meaning this",format(dataStartDate,"%B %Y"),

													ifelse(map_validation_df$county=="ka","Kauai County",

													       		 ifelse(map_validation_df$county=="oa","Honolulu County (Oahu)",

																       		 ifelse(map_validation_df$county=="mn","Maui County (Maui, Lanai, Molokai & Kahoolawe)",

																			       		 if(map_validation_df$county=="bi"){"Hawaii county"}))), 

							 		"monthly rainfall (mm) map is a", 

													ifelse(map_validation_df$rsq_rf_mm>=0.8,"high quality estimate of monthly rainfall.",

													       		 ifelse(map_validation_df$rsq_rf_mm>=0.6,"good quality estimate of monthly rainfall.",

																       		 ifelse(map_validation_df$rsq_rf_mm>=0.4,"moderate quality estimate of monthly rainfall.",

																			       		 ifelse(map_validation_df$rsq_rf_mm>=0.2,"low quality estimate of monthly rainfall, and should be used with dilligence.","lowest quality estimate of monthly rainfall, and should be used with caution.")))),

							 		"All maps are subject to change as new data becomes available or unknown errors are corrected in reoccurring versions.", 

									      	"Errors in rainfall estimates do vary over space meaning any gridded rainfall value, even on higher quality maps, could still produce incorrect estimates.",

											"Check standard error (SE) maps to better understand spatial estimates of prediction error" 

											))#statement end

					

				keywords<-as.character(paste(ifelse(map_validation_df$county=="ka","Kauai,",

								    			ifelse(map_validation_df$county=="oa","Oahu,",

												   		 	ifelse(map_validation_df$county=="mn","Maui, Lanai, Molokai, Kahoolawe,)",

																   		 	if(map_validation_df$county=="bi"){"Hawaii Island,"}))),"Hawaii, Rainfall prediction, Monthly precipitation, Rainfall, Climate, Spatial interpolation, Kriging"))





				file_meta<-data.frame(dataStatement=co_statement,

						      		keywords=keywords,

										county=map_validation_df$county,

										dataYearMon=gsub("_","-",map_validation_df$moYr),

												dataStartDate=dataStartDate,

												dataEndDate=dataEndDate,

														dateProduced=Sys.Date(),

														dataVersionType=datatype,

																RFstationFile=filenames[1],

																RFmmGridFile=filenames[2],

																		RFmmSEgridFile=filenames[3],

																		RFanomGridFile=filenames[4],

																				RFanomSEgridFile=filenames[5],

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

												   				value=as.character(c("All data produced by University of Hawaii at Manoa Dept. of Geography and the Enviroment, Ecohydology Lab in collaboration with the Water Resource Research Center (WRRC). Support for the Hawai?i EPSCoR Program is provided by the National Science Foundation?s Research Infrastructure Improvement (RII) Track-1: ?Ike Wai: Securing Hawaii?s Water Future Award # OIA-1557349.",

																						     "Matthew Lucas (mplucas@hawaii.edu), Ryan Longman (rlongman@hawaii.edu),Thomas Giambelluca (thomas@hawaii.edu)")))		

											)

					final_meta$attribute<-as.character(final_meta$attribute)

						final_meta<-final_meta[-c(grep("Class",final_meta$attribute),grep("probUgly",final_meta$attribute),grep("level",final_meta$attribute)),]#remove attribute rows

						row.names(final_meta)<-NULL

							final_meta[c(30:37),1]<-c("rsqRFmm","rmseRFmm","maeRFmm","biasRFmm","rsqRFanom","rmseRFanom","maeRFanom","biasRFanom")

							return(final_meta)

								}# county meta



			if(statewide){

							rf_validation<-data.frame(

										  				stationCountCounties=paste(map_validation_df$stationCount,collapse=", "),

																		krigeRFmmMaxCounties=paste(round(map_validation_df$krigeRFmmMax,5),collapse=", "),

																		krigeRFmmMinCounties=paste(round(map_validation_df$krigeRFmmMin,5),collapse=", "),

																						negRFgridZeroCounties=paste(map_validation_df$negRFgridZero,collapse=", "),

																						negRFpixCountCounties=paste(map_validation_df$negRFpixCount,collapse=", "),

																										gridPixCountCounties=paste(map_validation_df$gridPixCount,collapse=", "),

																										negRFloocvCounties=paste(map_validation_df$negRFloocv,collapse=", "),

																														rsqRFmmCounties=paste(round(map_validation_df$rsq_rf_mm,5),collapse=", "),

																														rmseRFmmCounties=paste(round(map_validation_df$rmse_rf_mm,5),collapse=", "),

																																		maeRFmmCounties=paste(round(map_validation_df$mae_rf_mm,5),collapse=", "),

																																		biasRFmmCounties=paste(round(map_validation_df$bias_rf_mm,5),collapse=", "),

																																						rsqRFanomCounties=paste(round(map_validation_df$rsq_rf_anom,5),collapse=", "),

																																						rmseRFanomCounties=paste(round(map_validation_df$rmse_rf_anom,5),collapse=", "),

																																										maeRFanomCounties=paste(round(map_validation_df$mae_rf_anom,5),collapse=", "),

																																										biasRFanomCounties=paste(round(map_validation_df$bias_rf_anom,5),collapse=", "),

																																														stationCount=sum(map_validation_df$stationCount),

																																														krigeRFmmMax=round(maxValue(hi_statewide_rf_mm_ras),5),

																																																		krigeRFmmMin=round(min(map_validation_df$krigeRFmmMin),5),

																																																		negRFgridZero=as.logical(max(map_validation_df$negRFgridZero)),

																																																						negRFpixCount=sum(map_validation_df$negRFpixCount),

																																																						gridPixCount=sum(map_validation_df$gridPixCount),

																																																										negRFloocv=as.logical(max(map_validation_df$negRFloocv)),

																																																										rsqRFmm=round(summary(lm(loocv_df$obs_rf~loocv_df$pred_rf))$r.squared,5),

																																																														rmseRFmm=round(rmse(loocv_df$obs_rf,loocv_df$pred_rf),5),

																																																														maeRFmm=round(mae(loocv_df$obs_rf,loocv_df$pred_rf),5),

																																																																		biasRFmm=round(bias(loocv_df$obs_rf,loocv_df$pred_rf),5),

																																																																		rsqRFanom=round(summary(lm(loocv_df$obs_anom~loocv_df$pred_anom))$r.squared,5),

																																																																						rmseRFanom=round(rmse(loocv_df$obs_anom,loocv_df$pred_anom),5),

																																																																						maeRFanom=round(mae(loocv_df$obs_anom,loocv_df$pred_anom),5),

																																																																										biasRFanom =round(bias(loocv_df$obs_anom,loocv_df$pred_anom),5)

																																																																										)





				state_statement<-as.character(paste("This",format(dataStartDate,"%B %Y"),"mosaic rainfall map of the State of Hawaii",

								    		"is a high spatial resolution (~250m) gridded prediction of cumulative rainfall in millimeters from",

												format(dataStartDate,"%d %b %Y"),"to",format(dataEndDate,"%d %b %Y."),

												"This was produced using several climate-aided modified automatic kriging interpolations of a log transformed sum of a monthly rainfall anomaly ratio (observed mm / mean monthly mm) and a constant of 1 (k), for each county extent. This process was done for four individually produced maps of Kauai, Honolulu (Oahu), Maui (Maui, Lanai, Molokai, & Kahoolawe) and Hawaii counties.", 

														"These kriging processes used", 

														sum(map_validation_df$stationCount), "unique station locations statewide",

																"and their",format(dataStartDate,"%B %Y"), 

																"recorded and/or estimated rainfall (mm) totals. Please consult each county map meta-data files for more details on map production and accuracy at the county scale.",  

																		"A leave one out cross validation (LOOCV) of the all station data used in all four counties produced individual R-squared values of:",

																		paste(round(map_validation_df[map_validation_df$county=="ka","rsq_rf_mm"],2),

																		      			round(map_validation_df[map_validation_df$county=="oa","rsq_rf_mm"],2),

																								round(map_validation_df[map_validation_df$county=="mn","rsq_rf_mm"],2),

																								round(map_validation_df[map_validation_df$county=="bi","rsq_rf_mm"],2),collapse=", "),

								    		"for Kauai, Honolulu (Oahu), Maui (Maui, Lanai, Molokai, & Kahoolawe) and Hawaii counties respectively.",

												"As a whole leave one out cross validation (LOOCV) data from all stations data compared to observed monthly rainfall (mm) produces a statewide R-squared value of:", 

												round(rf_validation$rsqRFmm,2), "meaning overall this",format(dataStartDate,"%B %Y"),"statewide mosaic monthly rainfall (mm) map is a",  

														ifelse(map_validation_df$rsq_rf_mm>=0.8,"high quality estimate of monthly rainfall.",

														       		 ifelse(map_validation_df$rsq_rf_mm>=0.6,"good quality estimate of monthly rainfall.",

																	       		 ifelse(map_validation_df$rsq_rf_mm>=0.4,"moderate quality estimate of monthly rainfall.",

																				       		 ifelse(map_validation_df$rsq_rf_mm>=0.2,"low quality estimate of monthly rainfall, and should be used with dilligence.","lowest quality estimate of monthly rainfall, and should be used with caution.")))),

								    		"All maps are subject to change as new data becomes available or unknown errors are corrected in reoccurring versions.", 

										      	"Errors in rainfall estimates do vary over space meaning any gridded rainfall value, even on higher quality maps, could still produce incorrect estimates.",

												"Check standard error (SE) maps to better understand spatial estimates of prediction error." 

												))#statewide statement end

						

					keywords<-"Hawaii, Hawaiian Islands, Rainfall prediction, Monthly precipitation, Rainfall, Climate, Spatial interpolation, Kriging"

					

				file_meta<-data.frame(dataStatement=state_statement,

						      		keywords=keywords,

										county=paste(map_validation_df$county,collapse=", "),

										dataYearMon=gsub("_","-",map_validation_df$moYr[1]),

												dataStartDate=dataStartDate,

												dataEndDate=dataEndDate,

														dateProduced=Sys.Date(),

														dataVersionType=datatype,

																RFstationFile=filenames[1],

																RFmmGridFile=filenames[2],

																		RFmmSEgridFile=filenames[3],

																		RFanomGridFile=filenames[4],

																				RFanomSEgridFile=filenames[5],

																				GeoCoordUnits="Decimal Degrees",

																						GeoCoordRefSystem=as.character(crs(grid)),

																						Xresolution=xres(grid),

																								Yresolution=yres(grid),

																								ExtentXmin=extent(grid)[1],

																										ExtentXmax=extent(grid)[2],

																										ExtentYmin=extent(grid)[3],

																												ExtentYmax=extent(grid)[4],

																												varioCounties=paste(map_validation_df$vario,collapse=", ")

																														)



					final_meta<-rbind(data.frame(attribute=as.character(names(file_meta)),value=as.character(lapply(file_meta[1,], dataconvert))),

							  				data.frame(attribute=as.character(names(rf_validation)),value=as.character(lapply(rf_validation[1,],dataconvert))),

															data.frame(attribute=as.character(c("credits", "contacts")),

																   				value=as.character(c("All data produced by University of Hawaii at Manoa Dept. of Geography and the Enviroment, Ecohydology Lab in collaboration with the Water Resource Research Center (WRRC). Support for the Hawai?i EPSCoR Program is provided by the National Science Foundation?s Research Infrastructure Improvement (RII) Track-1: ?Ike Wai: Securing Hawaii?s Water Future Award # OIA-1557349.",

																										     		"Matthew Lucas (mplucas@hawaii.edu), Ryan Longman (rlongman@hawaii.edu),Thomas Giambelluca (thomas@hawaii.edu)")))		

															)

					final_meta$attribute<-as.character(final_meta$attribute)

						return(final_meta)

						}#end statewide meta

			  }#end metamaker func


