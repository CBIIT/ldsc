FROM public.ecr.aws/amazonlinux/amazonlinux:2023

# Install necessary packages
RUN dnf -y update && \
    dnf -y install \
    gcc \
    gzip \
    mariadb105 \
    mariadb105-devel \
    pkgconf \
    pkgconf-pkg-config \
    bzip2 \
    wget \
    && dnf clean all

# Install Miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /miniconda.sh && \
    bash /miniconda.sh -b -p /opt/conda && \
    rm /miniconda.sh

# Add conda to PATH
ENV PATH=/opt/conda/bin:$PATH

# Create application directory
RUN mkdir /app

# Set the working directory
WORKDIR /app

# Copy environment.yml to the container
COPY environment3.yml /app/environment.yml

# Create the conda environment
RUN conda env create -f /app/environment.yml

# Ensure conda is initialized
RUN /opt/conda/bin/conda init bash

# Copy the rest of your application code to the container
COPY . /app

# Command to run your application with the environment activated
CMD ["bash", "-c", "source /opt/conda/etc/profile.d/conda.sh && conda activate ldscTest && exec bash"]


