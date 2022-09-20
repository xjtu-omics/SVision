#FROM python:3.6
#
#WORKDIR .
#
#ADD . .
#
#RUN pip install -r requirements.txt
#
#CMD ["python", "SVision"]

#FROM continuumio/miniconda3

#WORKDIR .
#
#COPY . .
#
#RUN conda env create -f environment.yml
#
## Make RUN commands use the new environment:
#SHELL ["conda", "run", "-n", "svisionenv", "/bin/bash", "-c"]
#
## Demonstrate the environment is activated:
#RUN echo "Make sure tensorflow is installed:"
#RUN python -c "import pyvcf"
#
#ENV PATH /opt/conda/envs/svisionenv/bin:$PATH

# The code to run when container is started:
#COPY SVision .
#ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "svisionenv", "python"]

#WORKDIR .
#COPY environment.yml .
#
#RUN conda env create -f environment.yml
#
#SHELL ["conda", "run", "-n", "svisionenv", "/bin/bash", "-c"]
#
#RUN conda init bash
#RUN echo "conda activate svisionenv" > ~/.bashrc
#
#RUN pip install SVision

#RUN conda init bash \
#    && conda activate sivsionenv \

## Install required packages to conda environment
#RUN pip install numpy==1.16.4 \
#    && pip install scipy==1.5.4 \
#    && pip install opencv-python-headless \
#    && pip install intervaltree \
#    && pip install beautifulsoup4 \
#    && pip install pysam \
#    && pip install pyvcf \
#    && pip install tensorflow==1.14.0 \
#    && pip install SVision
#CMD ["/opt/conda/envs/svisionenv/bin/SVision", "--help"]


FROM ubuntu:18.04
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install --assume-yes apt-utils
RUN apt-get update && apt-get install -y git python3>=3.6 python3-pip

RUN cd /tmp \
    && git clone https://github.com/xjtu-omics/SVision.git \
    && cd SVision \
    && pip3 install . \
    && cp -r /usr/local/bin/SVision /usr/bin/

