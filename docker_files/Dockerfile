# Get python docker image version 3.7
FROM python:3.7

# Create working directory variable
ENV WORKDIR=/submissions-pipeline

# Set up volume directory in docker
VOLUME ${WORKDIR}

# Set up working directory in docker
WORKDIR ${WORKDIR}

# Define a system argument
ARG DEBIAN_FRONTEND=noninteractive

# Install system libraries of general use and python v3.5 and odbc connector
RUN apt-get update && apt-get install -y --no-install-recommends build-essential \
    dos2unix

# Copy python requirements file to docker images
COPY requirements.txt /opt/public-repository-submissions/requirements.txt

# Install python requirements
RUN pip install --no-cache-dir -r /opt/public-repository-submissions/requirements.txt

# Make the app available at port 3838
EXPOSE 3838

# Make a data and test_input directory to store metadata and fasta files
RUN mkdir /opt/public-repository-submissions/data
RUN mkdir /opt/public-repository-submissions/test_input

# Copy files to docker images
COPY config_files/required_columns.yaml /opt/public-repository-submissions/config_files/required_columns.yaml
COPY config_files/default_config.yaml /opt/public-repository-submissions/config_files/default_config.yaml
COPY test_input/test_metadata.tsv /opt/public-repository-submissions/test_input/test_metadata.tsv

# Copy python scripts to docker images
COPY docker_files/extract_metadata.py /opt/public-repository-submissions/extract_metadata.py
COPY docker_files/extract_config.py /opt/public-repository-submissions/extract_config.py
COPY docker_files/extract_upload_log.py /opt/public-repository-submissions/extract_upload_log.py

COPY docker_files/update_config.py /opt/public-repository-submissions/update_config.py
COPY docker_files/update_upload_log.py /opt/public-repository-submissions/update_upload_log.py

COPY docker_files/check.py /opt/public-repository-submissions/check.py
COPY docker_files/retrieve.py /opt/public-repository-submissions/retrieve.py

# Copy pipeline scripts to docker images
COPY biosample_sra_submission.py /opt/public-repository-submissions/biosample_sra_submission.py
COPY genbank_submission.py /opt/public-repository-submissions/genbank_submission.py
COPY gisaid_submission.py /opt/public-repository-submissions/gisaid_submission.py
COPY gisaid_uploader.py /opt/public-repository-submissions/gisaid_uploader.py
COPY seqsender.py /opt/public-repository-submissions/seqsender.py
COPY submission_preparation.py /opt/public-repository-submissions/submission_preparation.py

# Copy bash scripts to docker images
COPY save_config.sh /opt/public-repository-submissions/save_config.sh
COPY public-repository-submissions.sh /opt/public-repository-submissions/public-repository-submissions.sh

# Convert public-repository-submissions.sh from Windows style line endings to Unix-like control characters
RUN dos2unix /opt/public-repository-submissions/public-repository-submissions.sh

# Allow permission to excute the bash scripts
RUN chmod a+x /opt/public-repository-submissions/public-repository-submissions.sh

# Allow permission to read, write, and execute files in submissions-pipeline directory
RUN chmod -R a+rwx ${WORKDIR}

# Clean up
RUN apt-get autoremove -y && apt-get clean -y && rm -rf /var/lib/apt/lists/*

# Execute the pipeline 
ENTRYPOINT ["/bin/bash", "/opt/public-repository-submissions/public-repository-submissions.sh"]








