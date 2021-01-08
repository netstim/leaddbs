#!/bin/bash
# My example bash script
#DIRECTORY=`dirname $0`
#echo $DIRECTORY

docker run --name OSS_container --volume $2:/opt/OSS-DBS --volume $1:/opt/Patient -it --rm custom_oss_platform python3 Launcher_OSS_lite.py
