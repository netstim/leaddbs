#!/bin/sh

# launch MATLAB
if [ "$ARCH" == "Linux" ]; then
    $ARTENOLIS_SOFT_PATH/MATLAB/$MATLAB_VER/bin/./matlab -nodesktop -nosplash < test/testAll.m

elif [ "$ARCH" == "macOS" ]; then
    caffeinate -u &
    /Applications/MATLAB_$MATLAB_VER.app/bin/matlab -nodisplay -nosplash < test/testAll.m

elif [ "$ARCH" == "windows" ]; then
    # change to the build directory
    echo " -- changing to the build directory --"
    cd "$ARTENOLIS_DATA_PATH\jenkins\\workspace\\$CI_PROJECT_NAME\\MATLAB_VER\\$MATLAB_VER\\label\\$NODE_LABELS"

    echo " -- launching MATLAB --"
    unset Path
    nohup "$ARTENOLIS_SOFT_PATH\MATLAB\\$MATLAB_VER\\\bin\\matlab.exe" -nojvm -nodesktop -nosplash -useStartupFolderPref -logfile output.log -wait -r "restoredefaultpath; cd $ARTENOLIS_DATA_PATH\jenkins\\workspace\\$CI_PROJECT_NAME\\MATLAB_VER\\$MATLAB_VER\\label\\$NODE_LABELS; cd test; testAll;" & PID=$!

    # follow the log file
    tail -n0 -F --pid=$! output.log 2>/dev/null

    # wait until the background process is done
    wait $PID
fi

CODE=$?
exit $CODE
