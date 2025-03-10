#
# to build this docker container image
# step 1: git clone https://github.com/GPlates/gplately.git
# step 2: cd gplately
# step 3: docker build -f docker/Dockerfile -t gplately .
#
# run the docker container
# docker run --rm -ti -p 8888:8888 -v `pwd`:/workspace/my_stuff gplately 
#

# if error, try with --no-cache option

FROM mambaorg/micromamba:1.5.8

LABEL org.opencontainers.image.authors="michael.chin@sydney.edu.au"
LABEL version="v2.0.0"

USER root
RUN apt-get update && apt-get install ffmpeg libsm6 libxext6  -y

USER $MAMBA_USER
COPY --chown=$MAMBA_USER:$MAMBA_USER ./docker/env.yaml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1 

COPY --chown=$MAMBA_USER:$MAMBA_USER . /tmp/gplately
RUN pip3 install /tmp/gplately
RUN rm -rf /tmp/gplately

COPY --chown=$MAMBA_USER:$MAMBA_USER ./Notebooks /workspace/Notebooks

WORKDIR /workspace

EXPOSE 8888

CMD ["jupyter", "notebook", "--no-browser", "--ip=0.0.0.0", "--NotebookApp.token=''"]
