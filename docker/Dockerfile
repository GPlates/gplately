#
# to build this docker container image
# step 1: git clone https://github.com/GPlates/gplately.git
# step 2: cd gplately
# step 3: `docker build -f docker/Dockerfile -t gplately .`
#
# run the docker container
#   `docker run --rm -ti --memory=24g -p 8888:8888 -v `pwd`:/workspace/my_stuff gplately`
#

# If error, try with --no-cache option

# If the Jupyter Notebook kernel dies when running with Docker, it could be due to insufficient memory. Try reducing the number of parallel processes. Or you can increase the memory allocated to Docker with `--memory=24g` when starting the Docker container.

# You may also need to increase the memory allocated to Docker Desktop on your host machine. For example, on
# macOS (Docker Desktop):

# - Open Docker Desktop → Settings (gear icon) → Resources → Advanced.
# - Drag the Memory slider up (or type a value) — set it to something like 24G.
# - Click Apply & Restart 
# - Re-check: docker info | grep -i "Total Memory"

FROM mambaorg/micromamba:2.9.0

LABEL org.opencontainers.image.authors="michael.chin@sydney.edu.au"
LABEL version="v2.1.0"

USER root
RUN apt-get update && apt-get install ffmpeg libsm6 libxext6  -y

USER $MAMBA_USER
COPY --chown=$MAMBA_USER:$MAMBA_USER ./docker/env.yaml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1 

RUN git clone --depth=1 https://github.com/GPlates/gplately.git /tmp/gplately
#COPY --chown=$MAMBA_USER:$MAMBA_USER . /tmp/gplately
RUN pip3 install /tmp/gplately
RUN rm -rf /tmp/gplately

COPY --chown=$MAMBA_USER:$MAMBA_USER ./Notebooks /workspace/Notebooks

WORKDIR /workspace

EXPOSE 8888

CMD ["jupyter", "notebook", "--no-browser", "--ip=0.0.0.0", "--NotebookApp.token=''"]
