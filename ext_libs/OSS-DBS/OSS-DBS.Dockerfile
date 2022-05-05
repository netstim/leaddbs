FROM ubuntu:20.04
MAINTAINER "Ningfei Li" <ningfei.li@gmail.com>

ENV DEBIAN_FRONTEND=noninteractive

WORKDIR /tmp

# Install packages from Ubuntu source
RUN apt-get update && apt-get dist-upgrade -y && \
    apt-get install -y --no-install-recommends software-properties-common && \
    apt-add-repository ppa:fenics-packages/fenics && \
    apt-get install -y --no-install-recommends fenics libfreeimage3 libgl1-mesa-glx \
    libglu1-mesa libmeschach-dev libreadline-dev libtbb2 libxcursor1 make neuron-dev \
    python3-dev python3-distutils python3-h5py python3-matplotlib python3-neuron python3-nibabel \
    python3-pandas python3-progress python3-psutil python3-tk wget xauth xterm

# Install Gmsh (package provided in Ubuntu source is buggy)
RUN wget -c https://gmsh.info/bin/Linux/gmsh-4.10.1-Linux64.tgz 2> /dev/null && \
    tar -xzf gmsh-4.10.1-Linux64.tgz && cp -r gmsh-4.10.1-Linux64/* /usr/local

# Install SALOME
RUN wget -c https://files.salome-platform.org/Salome/Salome9.8.0/SALOME-9.8.0-native-UB20.04-SRC.tar.gz 2> /dev/null && \
    tar -xzf SALOME-9.8.0-native-UB20.04-SRC.tar.gz && mv /tmp/SALOME-9.8.0-native-UB20.04-SRC /usr/local/salome && \
    sed -i 's/arch-linux-c-opt//' /usr/local/salome/salome && \
    rm -rf /usr/local/salome/SOURCES && rm -rf /usr/local/salome/ARCHIVES && rm -rf /tmp/*

# Update PATH
ENV PATH=/usr/local/salome:/usr/lib/neuron/bin:$PATH

# Suppress OpenMPI warning
RUN echo "btl_base_warn_component_unused = 0" >> /etc/openmpi/openmpi-mca-params.conf

# Set time zone
ENV TZ=Europe/Berlin

# Set default working directory
WORKDIR /opt/OSS-DBS/OSS_platform

CMD [ "/bin/bash" ]
