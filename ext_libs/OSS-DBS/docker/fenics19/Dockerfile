# Dockerfile to build the OSS DBS environment
#
# Author:
# Max Schroeder <max.schroeder@uni-rostock.de>


###############################################################################
# Stage 1
FROM sfbelaine/oss_dbs:base_latest as fenics_dep

WORKDIR /tmp

# Update current instance
RUN apt update && \
    apt full-upgrade -y && \
    apt clean

# Add fenics PPA repo for latest version
RUN apt update && \
    apt install -y \
        software-properties-common && \
    apt clean && \
    rm -rf /var/lib/apt/lists/*
RUN add-apt-repository ppa:fenics-packages/fenics

# Install dependencies available via apt-get
RUN apt update && \
    apt install -y \
        fenics \
        python3-dolfin \
        libgl1-mesa-glx \
        libxi-dev \
        libxmu-dev \
        libglu1-mesa-dev \
        libxinerama1 \
        python3-nibabel \
        libxcursor-dev && \
    apt clean && \
    rm -rf /var/lib/apt/lists/*

CMD ["/bin/bash"]
