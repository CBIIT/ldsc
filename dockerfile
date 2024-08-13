# Use an appropriate base image
FROM continuumio/miniconda3

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

SHELL ["bash", "-c"]
RUN echo "source /opt/conda/etc/profile.d/conda.sh && conda activate ldscTest" >> ~/.bashrc

CMD ["bash"]
# Command to run your application with the environment activated
#CMD ["bash", "-c", "source /opt/conda/etc/profile.d/conda.sh && conda activate ldscTest && exec bash"]