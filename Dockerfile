
############# Build Stage: Dependencies ##################

# Start from a base image
FROM --platform=linux/amd64 ubuntu:focal as base

# Install system libraries of general use
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update --allow-releaseinfo-change \
    && apt-get install --no-install-recommends -y \
    build-essential \ 
    python3.8 \
    python3-pip \
    awscli \
    groff \
    dpkg \
    sudo \
    curl \
    unzip \
    xtail \
    dos2unix \
    && apt clean autoclean \
    && apt autoremove --yes \
    && rm -rf /var/lib/{apt,dpkg,cache,log}/

############# Build Stage: Final ##################

# Build the final image 
FROM base as final

# Create program directory 
ENV PROGRAM=/seqsender

# Create working directory 
ENV WORKDIR=/data

# Set up volume directory in docker
VOLUME ${WORKDIR}

# Set up working directory in docker
WORKDIR ${WORKDIR}

############# Program Scripts ##################

# Copy all scripts to docker images
COPY . ${PROGRAM}

############# Install AWS CLI ##################

RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "${PROGRAM}/aws_cli/awscliv2.zip" && \
    unzip ${PROGRAM}/aws_cli/awscliv2.zip && \
    sudo ./aws/install    
    
RUN rm -rf ${PROGRAM}/aws_cli

RUN mkdir /root/.aws

RUN chmod 755 /root/.aws

ENV PATH=/usr/local/bin/aws:"${PATH}"

############# Install GISAID EpiFlu CLI ##################

# Copy all files to docker images
COPY epiflu_cli ${PROGRAM}

############# Install GISAID EpiCov CLI ##################

# Copy all files to docker images
COPY epicov_cli ${PROGRAM}

############# Install python packages ##################

# Copy python requirements file to docker images
COPY requirements.txt ${PROGRAM}/requirements.txt

# Install python requirements
RUN pip3 install --no-cache-dir -r ${PROGRAM}/requirements.txt

############# Launch program ##################

# Copy bash script to start the program to docker images
COPY seqsender-kickoff ${PROGRAM}/seqsender-kickoff

# Convert bash script from Windows style line endings to Unix-like control characters
RUN dos2unix ${PROGRAM}/seqsender-kickoff

# Allow permission to excute the bash scripts
RUN chmod a+x ${PROGRAM}/seqsender-kickoff

# Export bash script to path
ENV PATH "$PATH:${PROGRAM}"
