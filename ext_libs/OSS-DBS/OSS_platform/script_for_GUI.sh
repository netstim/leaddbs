#!/bin/bash
# My example bash script
#DIRECTORY=`dirname $0`
#echo $DIRECTORY
#python3 GUI_tree_files/AppUI.py $1

#Terminal.app -e python3 GUI_tree_files/AppUI.py $1
#docker run --volume $HOME/OSS-DBS:/opt/OSS-DBS --volume --cap-add=SYS_PTRACE -it --rm custom_oss_platform python3 Launcher_OSS_lite.py
#docker run --volume $HOME/OSS-DBS:/opt/OSS-DBS --volume $1 -it --rm sfbelaine/oss_dbs:platform_latest python3 Launcher_OSS_lite.py
#var1=$1
#var2=$2
#echo $var1
#echo $var2
#osascript -e 'tell application "Terminal" to activate' -e 'tell application "Terminal" to do script "cd Documents/MATLAB_files/leaddbs-oss-dbs/ext_libs/OSS-DBS/OSS_platform/"' -e 'tell application "Terminal" to do script "python3 GUI_tree_files/AppUI.py \"$1\""'
#osascript -e 'tell application "Terminal" to activate' -e 'tell application "Terminal" to do script "cd  \"$2\"" & " ; " & " python3 GUI_tree_files/AppUI.py \"$1\""'
osascript -e "tell application \"Terminal\" to do script \"cd $2; python3 GUI_tree_files/AppUI.py $1\""
