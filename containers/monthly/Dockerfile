# Start with the Air Temperature Preliminary Aggregation Base.
FROM ghcr.io/hcdp/task-preliminary-rainfall-base:latest
LABEL org.opencontainers.image.source="https://github.com/hcdp/preliminary_rainfall"
LABEL org.opencontainers.image.description="Monthly Air Temperature Aggregation and Mapping."

# Pull the monthly dependencies and extract them.
WORKDIR /home/hawaii_climate_products_container/preliminary/rainfall/dependencies/monthly
RUN wget https://ikeauth.its.hawaii.edu/files/v2/download/public/system/ikewai-annotated-data/HCDP/workflow_data/preliminary/rainfall/dependencies/monthly_dependencies.zip
RUN unzip monthly_dependencies.zip
RUN rm monthly_dependencies.zip

RUN mkdir -p /home/hawaii_climate_products_container/preliminary/rainfall/data_outputs
ADD /data_outputs /home/hawaii_climate_products_container/preliminary/rainfall/data_outputs

RUN mkdir -p /home/hawaii_climate_products_container/preliminary/rainfall/working_data
ADD /working_data /home/hawaii_climate_products_container/preliminary/rainfall/working_data

RUN mkdir -p /home/hawaii_climate_products_container/preliminary/rainfall/code/source
ADD /code/source /home/hawaii_climate_products_container/preliminary/rainfall/code/source

ADD /tools /sync


RUN mkdir -p /home/hawaii_climate_products_container/preliminary/rainfall/code/monthly
ADD /code/monthly /home/hawaii_climate_products_container/preliminary/rainfall/code/monthly

# Add the task wrapper script.
ADD /containers/monthly/task.sh /actor/task.sh

# Get the container's upload configs.
ADD /containers/monthly/configs/upload /sync

# This container is intended to run the task wrapper at launch, unless otherwise specified at run-time.
CMD [ "/bin/bash", "/actor/task.sh" ]