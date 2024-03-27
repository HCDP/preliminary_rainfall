###BASH CODE TO DEFINE & WGET ALL FILES IF THEY EXIST 

#define TZ as HST
TZ=Pacific/Honolulu 

#define filedate year_month of yesterday
if [ $AGGREGATION_DATE ]; then
    fileDate=$(date --date="$AGGREGATION_DATE" +"%Y_%m")
    fileYear=$(date --date="$AGGREGATION_DATE" +"%Y")
else
    fileDate=$(date --date="1 month ago" +"%Y_%m") 
    fileYear=$(date --date="1 month ago" +"%Y") 
fi

echo "fileDate is $fileDate"
echo "fileYear is $fileYear"

#define master dir
pathMaster=$'preliminary/rainfall/' 

#define url
urlMaster=$'https://ikeauth.its.hawaii.edu/files/v2/download/public/system/ikewai-annotated-data/HCDP/workflow_data/'

#define files
year_monthly_rf=$pathMaster'data_outputs/tables/station_data/monthly/partial_filled/statewide/Statewide_Partial_Filled_Monthly_RF_mm_'$fileYear'.csv'
month_daily_rf=$pathMaster'data_outputs/tables/station_data/daily/partial_filled/statewide/Statewide_Partial_Filled_Daily_RF_mm_'$fileDate'.csv'
month_count_rf=$pathMaster'data_outputs/tables/rf_station_tracking/count/monthly/'$fileYear'_count_log_monthly_rf.csv'

#set dir
cd /home/hawaii_climate_products_container

#wget files and put in proper dir locations
echo "Getting monthly data of this year from $urlMaster$year_monthly_rf"
echo "Getting daily data of this month from $urlMaster$month_daily_rf"
echo "Getting count from $urlMaster$month_count_rf"
wget $urlMaster$year_monthly_rf -O './'$year_monthly_rf || rm $year_monthly_rf
wget $urlMaster$month_daily_rf -O './'$month_daily_rf || rm $month_daily_rf
wget $urlMaster$month_count_rf -O './'$month_count_rf || rm $month_count_rf


##PAU##
