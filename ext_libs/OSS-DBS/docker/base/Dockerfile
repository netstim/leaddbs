# Dockerfile to build the OSS DBS environment derived from EMStim Tools:
# https://github.com/j-zimmermann/EMStimTools (GPLv3)
#
# This file uses the multi stage build feature of the docker service:
# https://docs.docker.com/develop/develop-images/multistage-build/
#
# Specifically the following stages are contained in that order extending the
# previous stage
# 1. base_dep
# 2. salome_dep
# 3. gmsh_dep
# 4. neuron_dep
# 5. paraview_dep
#
# Author:
# Max Schroeder <max.schroeder@uni-rostock.de>


###############################################################################
# Stage 1
FROM ubuntu:18.04 as base_dep

WORKDIR /tmp

# Update current instance
RUN apt update && \
    apt full-upgrade -y && \
    apt clean && \
    rm -rf /var/lib/apt/lists/*

# Install dependencies available via apt-get
RUN apt update && \
    apt install -y \
        software-properties-common \
        libgfortran3 \
        git \
        wget \
        bash-completion && \
    apt clean && \
    rm -rf /var/lib/apt/lists/*

# Install pip
RUN apt update && \
    apt install -y \
        python-pip \
        python3 \
        python3-pip && \
    apt clean && \
    rm -rf /var/lib/apt/lists/*

CMD ["/bin/bash"]


###############################################################################
# Stage 2
FROM base_dep as salome_dep

# Install SALOME
RUN wget -O SALOME-8.3.0-UB16.04.tgz 'https://www.salome-platform.org/downloads/previous-versions/salome-v8.3.0/DownloadDistr?platform=UB16.04&version=8.3.0' 2> /dev/null && \
    tar -xzf SALOME-8.3.0-UB16.04.tgz -C /opt && \
    rm -rf /tmp/*

# Update PATH
ENV PATH="${PATH}:/opt/SALOME-8.3.0-UB16.04"


###############################################################################
# Stage 3
FROM salome_dep as gmsh_dep

# Install GMSH and dependencies
RUN apt update && \
    apt install -y \
        libxft2 && \
    apt clean && \
    rm -rf /var/lib/apt/lists/*

RUN wget http://gmsh.info/bin/Linux/gmsh-4.2.2-Linux64-sdk.tgz 2> /dev/null && \
    tar -xzf gmsh-4.2.2-Linux64-sdk.tgz -C /opt && \
    rm -rf /tmp/*

# Update PATHs
ENV PATH="${PATH}:/opt/gmsh-4.2.2-Linux64-sdk/bin"
ENV PYTHONPATH="/opt/gmsh-4.2.2-Linux64-sdk/lib"


###############################################################################
# Stage 4
FROM gmsh_dep as neuron_dep

# Install NEURON
RUN apt update && \
    apt install -y \
        libx11-dev \
        bison \
        flex \
        automake \
        libtool \
        libxext-dev \
        libncurses-dev \
        python3-dev \
        xfonts-100dpi \
        cython3 \
        libopenmpi-dev \
        python3-scipy \
        make \
        zlib1g-dev && \
    apt clean && \
    rm -rf /var/lib/apt/lists/*
RUN wget https://neuron.yale.edu/ftp/neuron/versions/v7.6/iv-19.tar.gz 2> /dev/null && \
    tar -xzf iv-19.tar.gz -C /opt && \
    cd /opt/iv-19 && \
    ./configure && \
    make -j && \
    make install -j
RUN wget https://neuron.yale.edu/ftp/neuron/versions/v7.6/nrn-7.6.tar.gz 2> /dev/null && \
    tar -xzf nrn-7.6.tar.gz -C /opt && \
    cd /opt/nrn-7.6 && \
    ./configure --with-iv --with-paranrn --with-nrnpython=python3 && \
    make -j && \
    make install -j && \
    cd src/nrnpython && \
    python3 setup.py install && \
    rm -rf /tmp/*

# Update PATH
ENV PATH=$PATH:/usr/local/nrn/x86_64/bin/:/usr/local/iv/x86_64/bin/


###############################################################################
# Stage 5
FROM neuron_dep as paraview_dep

# Install ParaView
# Note that in order to call paraview from python scripts, at least the following environment variables have to be set:
# PYTHONPATH="/opt/ParaView-5.4.1-Qt5-OpenGL2-MPI-Linux-64bit/lib/python2.7/site-packages/:/opt/ParaView-5.4.1-Qt5-OpenGL2-MPI-Linux-64bit/lib/python2.7/site-packages/paraview/vtk/"
# LD_LIBRARY_PATH="/opt/ParaView-5.4.1-Qt5-OpenGL2-MPI-Linux-64bit/lib/paraview-5.4/"
RUN wget -O ParaView-5.4.1-Qt5-OpenGL2-MPI-Linux-64bit.tar.gz 'https://www.paraview.org/paraview-downloads/download.php?submit=Download&version=v5.4&type=binary&os=Linux&downloadFile=ParaView-5.4.1-Qt5-OpenGL2-MPI-Linux-64bit.tar.gz' 2> /dev/null && \
    tar -xzf ParaView-5.4.1-Qt5-OpenGL2-MPI-Linux-64bit.tar.gz -C /opt && \
    rm -rf /tmp/*

# Update PATH
ENV PATH="${PATH}:/opt/ParaView-5.4.1-Qt5-OpenGL2-MPI-Linux-64bit/bin"
