#!/bin/bash
# My example bash script
#DIRECTORY=`dirname $0`
#echo $DIRECTORY

docker run --volume $2:/opt/OSS-DBS --volume $1:/opt/Patient -it --rm ningfei/oss-dbs python3 Launcher_OSS_lite.py
