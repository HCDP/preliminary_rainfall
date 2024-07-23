#!/bin/bash
echo "[task.sh] [1/5] Starting Execution."
export TZ="HST"
echo "It is currently $(date)."
if [ $CUSTOM_DATE ]; then
    echo "An aggregation date was provided by the environment."
else
    export CUSTOM_DATE=$(date -d "1 day ago" --iso-8601)
    echo "No aggregation date was provided by the environment. Defaulting to yesterday."
fi
echo "Aggregation date is: " $CUSTOM_DATE

echo "[task.sh] [2/5] Acquiring yesterday's cumulative aggregation files, if exists."
echo "---begin daily_rf_wget.sh---"
bash /home/hawaii_climate_products_container/preliminary/rainfall/code/daily/daily_rf_wget.sh
echo "---end daily_rf_wget.sh---"

echo "[task.sh] [3/5] Aggregating Rainfall data on the daily timeframe."
echo "---begin hads_daily_rf_FINAL.R---"
Rscript /home/hawaii_climate_products_container/preliminary/rainfall/code/daily/hads_daily_rf_FINAL.R $CUSTOM_DATE
echo "---end hads_daily_rf_FINAL.R---"

echo "---begin nws_rr5_daily_rf_FINAL.R---"
Rscript /home/hawaii_climate_products_container/preliminary/rainfall/code/daily/nws_rr5_daily_rf_FINAL.R $CUSTOM_DATE
echo "---end nws_rr5_daily_rf_FINAL.R---"

echo "---begin madis_daily_rf_FINAL.R---"
Rscript /home/hawaii_climate_products_container/preliminary/rainfall/code/daily/madis_daily_rf_FINAL.R $CUSTOM_DATE
echo "---end madis_daily_rf_FINAL.R---"

echo "---begin himesoSyno_daily_rf_FINAL.R---"
Rscript /home/hawaii_climate_products_container/preliminary/rainfall/code/daily/himesoSyno_daily_rf_FINAL.R $CUSTOM_DATE
echo "---end himesoSyno_daily_rf_FINAL.R---"

echo "---begin all_data_daily_merge_table_rf_FINAL.R---"
Rscript /home/hawaii_climate_products_container/preliminary/rainfall/code/daily/all_data_daily_merge_table_rf_FINAL.R $CUSTOM_DATE
echo "---end all_data_daily_merge_table_rf_FINAL.R---"

echo "---begin qaqc_randfor_bad_data_flag_remove_rf_FINAL.R---"
Rscript /home/hawaii_climate_products_container/preliminary/rainfall/code/daily/qaqc_randfor_bad_data_flag_remove_rf_FINAL.R $CUSTOM_DATE
echo "---end qaqc_randfor_bad_data_flag_remove_rf_FINAL.R---"

echo "---begin daily_gap_fill_NR_only_rf_FINAL.R---"
Rscript /home/hawaii_climate_products_container/preliminary/rainfall/code/daily/daily_gap_fill_NR_only_rf_FINAL.R $CUSTOM_DATE
echo "---end daily_gap_fill_NR_only_rf_FINAL.R---"

echo "---begin all_data_daily_last_obs_FINAL.R---"
Rscript /home/hawaii_climate_products_container/preliminary/rainfall/code/daily/all_data_daily_last_obs_FINAL.R $CUSTOM_DATE
echo "---end all_data_daily_last_obs_FINAL.R---"

echo "---begin daily_rf_krig_map_makr_FINAL.R---"
Rscript /home/hawaii_climate_products_container/preliminary/rainfall/code/daily/daily_rf_krig_map_makr_FINAL.R $CUSTOM_DATE
echo "---end daily_rf_krig_map_makr_FINAL.R---"

echo "[task.sh] [4/5] Preparing upload config."
cd /sync
python3 inject_upload_config.py upload.json $CUSTOM_DATE

echo "[task.sh] [5/5] Uploading data."
python3 upload.py

echo "[task.sh] All done!"
