#!/bin/bash
###BASH CODE TO DEFINE & WGET ALL FILES IF THEY EXIST 

echo "Starting Daily RF Pull"
#define TZ as HST
TZ=Pacific/Honolulu 

#define fileDate year_month of yesterday
if [ $AGGREGATION_DATE ]; then
    fileDate=$(date --date="$AGGREGATION_DATE" +"%Y_%m")
else
    fileDate=$(date --date="1 day ago" +"%Y_%m")  
fi

#define master dir
pathMaster=$'preliminary/rainfall/' 

#define url
urlMaster=$'https://ikeauth.its.hawaii.edu/files/v2/download/public/system/ikewai-annotated-data/HCDP/workflow_data/'

#define files
NWS_daily_rf=$pathMaster'working_data/nws_rr5/'$fileDate'_nws_daily_rf.csv'
Hads_daily_rf=$pathMaster'working_data/hads/'$fileDate'_hads_daily_rf.csv'
Madis_daily_rf=$pathMaster'working_data/madis/'$fileDate'_madis_daily_rf.csv'
RF_source=$pathMaster'data_outputs/tables/station_data/daily/source/statewide/Statewide_Daily_Source_'$fileDate'.csv'
RF_raw=$pathMaster'data_outputs/tables/station_data/daily/raw/statewide/Statewide_Raw_Daily_RF_mm_'$fileDate'.csv'
QAQC_flag=$pathMaster'data_outputs/tables/station_data/daily/qc_flag/statewide/Statewide_RF_QAQC_Daily_Bad_Flag_'$fileDate'.csv'
QAQC_prob=$pathMaster'data_outputs/tables/station_data/daily/qc_prob/statewide/Statewide_RF_QAQC_Daily_Bad_Prob_'$fileDate'.csv'
QAQC_fail=$pathMaster'data_outputs/tables/rf_station_tracking/qaqc_fail/'$fileDate'_qaqc_fail_daily_rf.csv'
QAQC_raw=$pathMaster'data_outputs/tables/station_data/daily/raw_qc/statewide/Statewide_QAQC_Raw_Daily_RF_mm_'$fileDate'.csv'
RF_fill=$pathMaster'data_outputs/tables/station_data/daily/partial_filled/statewide/Statewide_Partial_Filled_Daily_RF_mm_'$fileDate'.csv'
Fill_log=$pathMaster'data_outputs/tables/rf_station_tracking/gapFilling/'$fileDate'_gap_fill_log_daily_rf.csv'
RF_unknown=$pathMaster'data_outputs/tables/rf_station_tracking/missing/'$fileDate'_unknown_rf_sta.csv'
Count_log=$pathMaster'data_outputs/tables/rf_station_tracking/count/daily/'$fileDate'_count_log_daily_rf.csv'
Last_RF=$pathMaster'data_outputs/tables/rf_station_tracking/lastObs/lastRFdayObs.csv'

#set dir
cd /home/hawaii_climate_products_container

#wget files and put in proper dir locations
wget $urlMaster$NWS_daily_rf -O './'$NWS_daily_rf || rm $NWS_daily_rf
wget $urlMaster$Hads_daily_rf -O './'$Hads_daily_rf || rm $Hads_daily_rf
wget $urlMaster$Madis_daily_rf -O './'$Madis_daily_rf || rm $Madis_daily_rf
wget $urlMaster$RF_source -O './'$RF_source || rm $RF_source
wget $urlMaster$RF_raw -O './'$RF_raw || rm $RF_raw
wget $urlMaster$QAQC_flag -O './'$QAQC_flag || rm $QAQC_flag   
wget $urlMaster$QAQC_prob -O './'$QAQC_prob || rm $QAQC_prob   
wget $urlMaster$QAQC_fail -O './'$QAQC_fail || rm $QAQC_fail   
wget $urlMaster$QAQC_raw -O './'$QAQC_raw || rm $QAQC_raw      
wget $urlMaster$RF_fill -O './'$RF_fill || rm $RF_fill
wget $urlMaster$Fill_log -O './'$Fill_log || rm $Fill_log      
wget $urlMaster$RF_unknown -O './'$RF_unknown || rm $RF_unknown
wget $urlMaster$Count_log -O './'$Count_log || rm $Count_log   
wget $urlMaster$Last_RF -O './'$Last_RF || rm $Last_RF

##PAU##