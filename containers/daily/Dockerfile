# Start with the Air Temperature Preliminary Aggregation Base.
FROM ghcr.io/hcdp/task-preliminary-rainfall-base:latest
LABEL org.opencontainers.image.source="https://github.com/hcdp/preliminary_rainfall"
LABEL org.opencontainers.image.description="Daily Air Temperature Aggregation and Mapping."

# Pull the daily dependencies and extract them.
RUN mkdir -p /home/hawaii_climate_products_container/preliminary/rainfall/dependencies
WORKDIR /home/hawaii_climate_products_container/preliminary/rainfall/dependencies
RUN wget -nv https://ikeauth.its.hawaii.edu/files/v2/download/public/system/ikewai-annotated-data/HCDP/workflow_data/preliminary/rainfall/dependencies/daily_dependencies.tar.gz
RUN tar -xf daily_dependencies.tar.gz
RUN rm daily_dependencies.tar.gz

RUN mkdir -p /home/hawaii_climate_products_container/preliminary/rainfall/data_outputs
ADD /data_outputs /home/hawaii_climate_products_container/preliminary/rainfall/data_outputs

RUN mkdir -p /home/hawaii_climate_products_container/preliminary/rainfall/working_data
ADD /working_data /home/hawaii_climate_products_container/preliminary/rainfall/working_data

RUN mkdir -p /home/hawaii_climate_products_container/preliminary/rainfall/code/source
ADD /code/source /home/hawaii_climate_products_container/preliminary/rainfall/code/source

ADD /tools /sync


RUN mkdir -p /home/hawaii_climate_products_container/preliminary/rainfall/code/daily
ADD /code/daily /home/hawaii_climate_products_container/preliminary/rainfall/code/daily

# Add upload config
ADD /containers/daily/configs/upload/upload_config.json /sync/upload_config.json

# Add the task wrapper script.
ADD /containers/daily/task.sh /actor/task.sh

# Get the container's upload configs.
ADD /containers/daily/configs/upload /sync

# This container is intended to run the task wrapper at launch, unless otherwise specified at run-time.
CMD [ "/bin/bash", "/actor/task.sh" ]