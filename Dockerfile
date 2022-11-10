
# Start from a base image
FROM --platform=linux/amd64 ubuntu:focal as base

# Define a system argument
ARG DEBIAN_FRONTEND=noninteractive

# Install system libraries of general use
RUN apt-get update --allow-releaseinfo-change && apt-get install --no-install-recommends -y \
    build-essential \ 
    python3.7\
    python3-pip \
    python3-setuptools \
    python3-dev \
    dos2unix

# Create working directory variable
ENV WORKDIR=/submissions-pipeline

# Set up volume directory in docker
VOLUME ${WORKDIR}

# Set up working directory in docker
WORKDIR ${WORKDIR}

# Allow permission to read and write files to current working directory
RUN chmod -R a+rwx ${WORKDIR}

############# Install python packages ##################

# Copy python requirements file to docker images
COPY requirements.txt /seqsender/requirements.txt

# Install python requirements
RUN pip3 install --no-cache-dir -r /seqsender/requirements.txt

############# Create template directory ##################

# Create a template directory to store sample config file and metadata worksheet
RUN mkdir /seqsender/template

# Copy template files to docker images
COPY config_files/default_config.yaml /seqsender/template/default_config_template.yaml
COPY config_files/required_columns.yaml /seqsender/template/required_columns_template.yaml
COPY test_input/test_metadata.tsv /seqsender/template/metadata_template.tsv

############# Create test_input directory ##################

# Create a test_input directory to store testing files
RUN mkdir /seqsender/test_input

# Copy submission.xml to test_input directory in order to run test_bioproject function
COPY test_input/submission.xml /seqsender/test_input/submission.xml

############# Create config_files directory ##################

# Create a config_files directory to store config files
RUN mkdir /seqsender/config_files

# Copy config files to docker images
COPY config_files/default_config.yaml /seqsender/config_files/default_config.yaml
COPY config_files/required_columns.yaml /seqsender/config_files/required_columns.yaml

############# Create data directory ##################

# Create a data directory to store metadata and fasta files
RUN mkdir /seqsender/data

############# Run seqsender ##################

# Copy all scripts to docker images
COPY . /seqsender

# Copy bash script files to docker images
COPY save_config.sh /seqsender/save_config.sh

# Convert save_config.sh from Windows style line endings to Unix-like control characters
RUN dos2unix /seqsender/save_config.sh

# Allow permission to excute the bash scripts
RUN chmod a+x /seqsender/save_config.sh

# Copy bash script files to docker images
COPY seqsender-kickoff /seqsender/seqsender-kickoff

# Convert save_config.sh from Windows style line endings to Unix-like control characters
RUN dos2unix /seqsender/seqsender-kickoff

# Allow permission to excute the bash scripts
RUN chmod a+x /seqsender/seqsender-kickoff

# Allow permission to read, write, and execute files in /seqsender directory
RUN chmod -R a+rwx /seqsender

# Clean up
RUN apt-get autoremove -y && apt-get clean -y && rm -rf /var/lib/apt/lists/*

# Export bash script to path
ENV PATH "$PATH:/seqsender"

# Execute the pipeline 
ENTRYPOINT ["bash", "seqsender-kickoff"]



