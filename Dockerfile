#FROM python:3.6
#
#WORKDIR .
#
#ADD . .
#
#RUN pip install -r requirements.txt
#
#CMD ["python", "SVision"]

FROM continuumio/miniconda3

WORKDIR .

COPY . .
RUN conda env create -f environment.yml

# Make RUN commands use the new environment:
SHELL ["conda", "run", "-n", "svisionenv", "/bin/bash", "-c"]

# Demonstrate the environment is activated:
RUN echo "Make sure tensorflow is installed:"
RUN python -c "import tensorflow"

ENV PATH /opt/conda/envs/svisionenv/bin:$PATH

# The code to run when container is started:
#COPY SVision .
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "svisionenv", "python"]