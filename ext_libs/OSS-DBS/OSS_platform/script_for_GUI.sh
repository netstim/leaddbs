#!/bin/bash
# My example bash script
#DIRECTORY=`dirname $0`
#echo $DIRECTORY

osascript -e "tell application \"Terminal\" to do script \"cd $2; python3 GUI_tree_files/AppUI.py $1 $3 $4\""
