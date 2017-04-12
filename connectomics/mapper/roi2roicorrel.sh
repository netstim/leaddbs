#!/bin/bash

##########################################################################################################
##  This program computes connectomic outputs using seed region based on standard connectomes           ##
##   script by Andreas Horn, August, 2016.                                                              ##
##                                                                                                      ##
##                                                                                                      ##
##########################################################################################################


usage() { # DEFINE usage function here; functions must be defined before they are invoked, so this subroutine needs to be initialized before we parse arguments below.

printf "\v\t This function runs connectivity between two rois"
}

# define defaults here
writesingle="0"
command="seed"
outputfolder=$HOME
maskname="."
doboth=1
dodMRI="0"
dofMRI="0"
dmriresolution="222"
fcname="GSP_1000_Yeo_2011"
dcname="Groupconnectome_Horn_2013"

[[ "$1" == "" ]] && usage && exit
while [ "$1" != "" ]; do
    case $1 in
        -h  | --help)           shift
                                usage #RUN USAGE HERE
                                exit
                                ;;
        -r | --roilist )        shift
                                roilistname=$1
                                ;;
        -o | --outfile )        shift
                                outfilename=$1
                                ;;
        -s | --sublist )        shift
                                sublistname=$1
                                ;;
    esac
    shift
done


        # setup command here:
        cmd="/autofs/cluster/nimlab/connectomes/software/lead_dbs/connectomics/mapper/run_roi2roi_correl.sh /usr/pubsw/common/matlab/8.6 $sublistname $roilistname $outfilename"

    # check if we are on launchpad:

    

        # run actual script:
        echo $cmd
        $cmd




function usage
{
    echo "usage: system_page [[[-r roi ] [-c command]] | [-h help]]"
}
