#!/bin/bash
echo "[task.sh] [1/7] Starting Execution."
export TZ="HST"
echo "It is currently $(date)."
if [ $AGGREGATION_DATE ]; then
    echo "An aggregation date was provided by the environment."
    echo "Aggregation date is: " $AGGREGATION_DATE
else
    export AGGREGATION_DATE=$(date -d "1 day ago" --iso-8601)
    echo "No aggregation date was provided by the environment. Defaulting to yesterday."
    echo "Aggregation date is: " $AGGREGATION_DATE
fi

echo "[task.sh] [2/7] Acquiring yesterday's cumulative aggregation files, if exists."
echo "---begin daily_rf_wget.sh---"
bash /home/hawaii_climate_products_container/preliminary/rainfall/code/daily/daily_rf_wget.sh
echo "---end daily_rf_wget.sh---"

echo "[task.sh] [3/7] Aggregating Rainfall data on the daily timeframe."
echo "---begin hads_daily_rf_FINAL.R---"
Rscript /home/hawaii_climate_products_container/preliminary/rainfall/code/daily/hads_daily_rf_FINAL.R $AGGREGATION_DATE
echo "---end hads_daily_rf_FINAL.R---"

echo "---begin nws_rr5_daily_rf_FINAL.R---"
Rscript /home/hawaii_climate_products_container/preliminary/rainfall/code/daily/nws_rr5_daily_rf_FINAL.R $AGGREGATION_DATE
echo "---end nws_rr5_daily_rf_FINAL.R---"

echo "---begin madis_daily_rf_FINAL.R---"
Rscript /home/hawaii_climate_products_container/preliminary/rainfall/code/daily/madis_daily_rf_FINAL.R $AGGREGATION_DATE
echo "---end madis_daily_rf_FINAL.R---"

echo "---begin himesoSyno_daily_rf_FINAL.R---"
Rscript /home/hawaii_climate_products_container/preliminary/rainfall/code/daily/himesoSyno_daily_rf_FINAL.R $AGGREGATION_DATE
echo "---end himesoSyno_daily_rf_FINAL.R---"

echo "---begin all_data_daily_merge_table_rf_FINAL.R---"
Rscript /home/hawaii_climate_products_container/preliminary/rainfall/code/daily/all_data_daily_merge_table_rf_FINAL.R $AGGREGATION_DATE
echo "---end all_data_daily_merge_table_rf_FINAL.R---"

echo "---begin qaqc_randfor_bad_data_flag_remove_rf_FINAL.R---"
Rscript /home/hawaii_climate_products_container/preliminary/rainfall/code/daily/qaqc_randfor_bad_data_flag_remove_rf_FINAL.R $AGGREGATION_DATE
echo "---end qaqc_randfor_bad_data_flag_remove_rf_FINAL.R---"

echo "---begin daily_gap_fill_NR_only_rf_FINAL.R---"
Rscript /home/hawaii_climate_products_container/preliminary/rainfall/code/daily/daily_gap_fill_NR_only_rf_FINAL.R $AGGREGATION_DATE
echo "---end daily_gap_fill_NR_only_rf_FINAL.R---"

echo "---begin all_data_daily_last_obs_FINAL.R---"
Rscript /home/hawaii_climate_products_container/preliminary/rainfall/code/daily/all_data_daily_last_obs_FINAL.R $AGGREGATION_DATE
echo "---end all_data_daily_last_obs_FINAL.R---"

echo "[task.sh] [4/7] Preparing for intermediate data upload."
cd /sync
python3 update_date_string_in_config.py intermediate_upload_config.json intermediate_upload_config_datestrings_loaded.json $AGGREGATION_DATE
python3 add_upload_list_to_config.py intermediate_upload_config_datestrings_loaded.json config.json
python3 add_auth_info_to_config.py config.json

echo "[task.sh] [5/7] Attempting to upload the aggregated intermediate data."
python3 upload.py

echo "[task.sh] [6/7] Preparing for production data upload."
cd /sync
python3 update_date_string_in_config.py upload_config.json upload_config_datestrings_loaded.json $AGGREGATION_DATE
python3 add_upload_list_to_config.py upload_config_datestrings_loaded.json config.json
python3 add_auth_info_to_config.py config.json

echo "[task.sh] [7/7] Attempting to upload the aggregated production data."
python3 upload.py

echo "[task.sh] All done!"