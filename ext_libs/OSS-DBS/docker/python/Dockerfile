# Dockerfile to build the OSS DBS environment
#
# Author:
# Max Schroeder <max.schroeder@uni-rostock.de>


###############################################################################
# Stage 1
FROM sfbelaine/oss_dbs:fenics19_latest as python_dep

WORKDIR /opt/OSS-DBS

# Install python packages
COPY requirements.txt /opt/OSS-DBS/requirements.txt
RUN pip3 install -r requirements.txt
