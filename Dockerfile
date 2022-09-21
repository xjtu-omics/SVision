FROM ubuntu:18.04
ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y cmake make build-essential liblzma-dev libbz2-dev zlib1g-dev libncurses5-dev libncursesw5-dev
RUN apt-get update && apt-get install --assume-yes apt-utils
RUN apt-get update && apt-get install -y git python3>=3.6 python3-pip python3-opencv


RUN git clone https://github.com/xjtu-omics/SVision.git \
    && cd SVision \
    && pip3 install scikit-build \
    && pip3 install Cython \
    && pip3 install protobuf==3.13.0 \
    && pip3 install . \
    && cp -r ./SVision /usr/bin/