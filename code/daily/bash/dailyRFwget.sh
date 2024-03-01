##wget dependencies files for daily rf run
cd /home/hawaii_climate_products_container/preliminary/rainfall/dependencies #set dir for dependencies
wget https://ikeauth.its.hawaii.edu/files/v2/download/public/system/ikewai-annotated-data/HCDP/rainfall/HCDP_dependicies/daily_dependencies.tar.gz #wget daily dependencies tar gz
tar -xf daily_dependencies.tar.gz  #extract dependencies dirs and files
rm daily_dependencies.tar.gz #remove dependencies tar gz

##wget datastream daily rf monthly files
wget https://ikeauth.its.hawaii.edu/files/v2/download/public/system/ikewai-annotated-data/HCDP/workflow_data/preliminary_test/rainfall/working_data/hads/"$(date --date="yesterday" +"%Y_%m")_hads_daily_rf.csv" -P /home/hawaii_climate_products_container/preliminary/rainfall/working_data/hads/
wget https://ikeauth.its.hawaii.edu/files/v2/download/public/system/ikewai-annotated-data/HCDP/workflow_data/preliminary_test/rainfall/working_data/madis/"$(date --date="yesterday" +"%Y_%m")_madis_daily_rf.csv" -P /home/hawaii_climate_products_container/preliminary/rainfall/working_data/madis/
wget https://ikeauth.its.hawaii.edu/files/v2/download/public/system/ikewai-annotated-data/HCDP/workflow_data/preliminary_test/rainfall/working_data/nws_rr5/"$(date --date="yesterday" +"%Y_%m")_nws_daily_rf.csv" -P /home/hawaii_climate_products_container/preliminary/rainfall/working_data/nws_rr5/

##wget station monitoring log tables
#missing/unknown stations monthly file
wget https://ikeauth.its.hawaii.edu/files/v2/download/public/system/ikewai-annotated-data/HCDP/workflow_data/preliminary_test/rainfall/data_outputs/tables/rf_station_tracking/missing/"$(date --date="yesterday" +"%Y_%m")_unknown_rf_sta.csv" -P /home/hawaii_climate_products_container/preliminary/rainfall/data_outputs/tables/rf_station_tracking/missing/
#count log 
wget https://ikeauth.its.hawaii.edu/files/v2/download/public/system/ikewai-annotated-data/HCDP/workflow_data/preliminary_test/rainfall/data_outputs/tables/rf_station_tracking/count/"$(date --date="yesterday" +"%Y_%m")_count_log_daily_rf.csv" -P /home/hawaii_climate_products_container/preliminary/rainfall/data_outputs/tables/rf_station_tracking/count/
#gap fill log
wget https://ikeauth.its.hawaii.edu/files/v2/download/public/system/ikewai-annotated-data/HCDP/workflow_data/preliminary_test/rainfall/data_outputs/tables/rf_station_tracking/gapFilling/"$(date --date="yesterday" +"%Y_%m")_gap_fill_log_daily_rf.csv" -P /home/hawaii_climate_products_container/preliminary/rainfall/data_outputs/tables/rf_station_tracking/gapFilling/
#qaqc fail log
wget https://ikeauth.its.hawaii.edu/files/v2/download/public/system/ikewai-annotated-data/HCDP/workflow_data/preliminary_test/rainfall/data_outputs/tables/rf_station_tracking/qaqc_fail/"$(date --date="yesterday" +"%Y_%m")_qaqc_fail_daily_rf.csv" -P /home/hawaii_climate_products_container/preliminary/rainfall/data_outputs/tables/rf_station_tracking/qaqc_fail/

##wget daily rf monthly station files 
#rf station daily source
wget https://ikeauth.its.hawaii.edu/files/v2/download/public/system/ikewai-annotated-data/HCDP/workflow_data/preliminary_test/rainfall/data_outputs/tables/station_data/daily/source/statewide/Statewide_Daily_Source_"$(date --date="yesterday" +"%Y_%m").csv" -P /home/hawaii_climate_products_container/preliminary/rainfall/data_outputs/tables/station_data/daily/source/statewide/
#rf station daily raw rainfall mm
wget https://ikeauth.its.hawaii.edu/files/v2/download/public/system/ikewai-annotated-data/HCDP/workflow_data/preliminary_test/rainfall/data_outputs/tables/station_data/daily/raw/statewide/Statewide_Raw_Daily_RF_mm_"$(date --date="yesterday" +"%Y_%m").csv" -P /home/hawaii_climate_products_container/preliminary/rainfall/data_outputs/tables/station_data/daily/raw/statewide/
#rf station daily qaqc raw rainfall mm
wget https://ikeauth.its.hawaii.edu/files/v2/download/public/system/ikewai-annotated-data/HCDP/workflow_data/preliminary_test/rainfall/data_outputs/tables/station_data/daily/raw_qc/statewide/Statewide_QAQC_Raw_Daily_RF_mm_"$(date --date="yesterday" +"%Y_%m").csv" -P /home/hawaii_climate_products_container/preliminary/rainfall/data_outputs/tables/station_data/daily/raw_qc/statewide/
#rf station daily qaqc flag
wget https://ikeauth.its.hawaii.edu/files/v2/download/public/system/ikewai-annotated-data/HCDP/workflow_data/preliminary_test/rainfall/data_outputs/tables/station_data/daily/qc_flag/statewide/Statewide_RF_QAQC_Daily_Bad_Flag_"$(date --date="yesterday" +"%Y_%m").csv" -P /home/hawaii_climate_products_container/preliminary/rainfall/data_outputs/tables/station_data/daily/qc_flag/statewide/
#rf station daily qaqc prob
wget https://ikeauth.its.hawaii.edu/files/v2/download/public/system/ikewai-annotated-data/HCDP/workflow_data/preliminary_test/rainfall/data_outputs/tables/station_data/daily/qc_prob/statewide/Statewide_RF_QAQC_Daily_Bad_Prob_"$(date --date="yesterday" +"%Y_%m").csv" -P /home/hawaii_climate_products_container/preliminary/rainfall/data_outputs/tables/station_data/daily/qc_prob/statewide/
#rf station daily partial filled rainfall mm
wget https://ikeauth.its.hawaii.edu/files/v2/download/public/system/ikewai-annotated-data/HCDP/workflow_data/preliminary_test/rainfall/data_outputs/tables/station_data/daily/partial_filled/statewide/Statewide_Partial_Filled_Daily_RF_mm_"$(date --date="yesterday" +"%Y_%m").csv" -P /home/hawaii_climate_products_container/preliminary/rainfall/data_outputs/tables/station_data/daily/partial_filled/


