#!/bin/bash
# My example bash script
#DIRECTORY=`dirname $0`
#echo $DIRECTORY
#python3 $1 I_am_with_stupid
#docker run --volume $HOME/OSS-DBS:/opt/OSS-DBS --volume --cap-add=SYS_PTRACE -it --rm custom_oss_platform python3 Launcher_OSS_lite.py
docker run --name OSS_container --volume $2:/opt/OSS-DBS --volume $1 -it --rm sfbelaine/oss_dbs:platform_latest python3 Launcher_OSS_lite.py
#osascript -e 'tell application "Terminal" to activate' -e 'tell application "Terminal" to do script "docker run --volume $HOME/OSS-DBS:/opt/OSS-DBS --volume $1 -it --rm sfbelaine/oss_dbs:platform_latest python3 $2/Launcher_OSS_lite.py"'
