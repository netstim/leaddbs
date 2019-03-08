#!/bin/sh

# launch MATLAB
if [ "$ARCH" == "Linux" ]; then
    $ARTENOLIS_SOFT_PATH/MATLAB/$MATLAB_VER/bin/./matlab -nodesktop -nosplash < test/testAll.m
fi

CODE=$?
exit $CODE
