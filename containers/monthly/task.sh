#!/bin/bash
echo "[task.sh] [1/7] Starting Execution."
export TZ="HST"
echo "It is currently $(date)."
if [ $CUSTOM_DATE ]; then
    export CUSTOM_DATE=$(date -d "$CUSTOM_DATE" +"%Y-%m-01")
    export CUSTOM_DATE=$(date -d "$CUSTOM_DATE +1 month -1 day" --iso-8601)
    echo "An aggregation date was provided by the environment. Setting to last day of provided month."
else
    export CUSTOM_DATE=$(date -d "-$(date +%d) days" --iso-8601)
    echo "No aggregation date was provided by the environment. Defaulting to last month."
fi

echo "Aggregation date is: " $CUSTOM_DATE

echo "[task.sh] [2/7] Acquiring Statewide Partially-filled Daily Rainfall data for this month."
cd /home/hawaii_climate_products_container/

echo "---begin monthly_rf_wget.sh---"
#uses CUSTOM_DATE env variable for date
bash /home/hawaii_climate_products_container/preliminary/rainfall/code/monthly/monthly_rf_wget.sh
echo "---end monthly_rf_wget.sh---"

echo "[task.sh] [3/7] Aggregating Rainfall data on the monthly timeframe."
cd /home/hawaii_climate_products_container/preliminary/rainfall/code/monthly

echo "---begin daily_to_monthly_agg_FINAL.R---"
Rscript /home/hawaii_climate_products_container/preliminary/rainfall/code/monthly/daily_to_monthly_agg_FINAL.R $CUSTOM_DATE
echo "---end daily_to_monthly_agg_FINAL.R---"

echo "---begin monthly_rf_krig_map_makr_FINAL.R---"
Rscript /home/hawaii_climate_products_container/preliminary/rainfall/code/monthly/monthly_rf_krig_map_makr_FINAL.R $CUSTOM_DATE
echo "---end monthly_rf_krig_map_makr_FINAL.R---"

echo "[task.sh] [4/7] Preparing to upload intermediate products."
cd /sync
python3 update_date_string_in_config.py intermediate_products.json intermediate_products_datestrings_loaded.json $CUSTOM_DATE
python3 add_upload_list_to_config.py intermediate_products_datestrings_loaded.json config.json
python3 add_auth_info_to_config.py config.json

echo "[task.sh] [5/7] Attempting to upload the intermediate files."
python3 upload.py
mv config.json intermediate_products_config.json

echo "[task.sh] [6/7] Preparing to upload final products."
cd /sync
python3 update_date_string_in_config.py final_products.json final_products_datestrings_loaded.json $CUSTOM_DATE
python3 add_upload_list_to_config.py final_products_datestrings_loaded.json config.json
python3 add_auth_info_to_config.py config.json

echo "[task.sh] [7/7] Attempting to upload the final files."
python3 upload.py

echo "[task.sh] All done!"