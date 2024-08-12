# syntax=docker/dockerfile:1

FROM debian:bullseye-slim

WORKDIR /app

RUN  apt-get update \
  && apt install -y bedtools python3-pip \
  && rm -rf /var/lib/apt/lists/* 

COPY . /app/ldsc/

RUN cd ldsc/ \
  && pip install .

CMD ["bash"]