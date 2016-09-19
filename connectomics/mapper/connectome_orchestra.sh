#!/bin/bash

##########################################################################################################
##  This program computes connectomic outputs using seed region based on standard connectomes           ##
##   script by Andreas Horn, August, 2016.                                                              ##
##                                                                                                      ##
##                                                                                                      ##
##########################################################################################################


usage() { # DEFINE usage function here; functions must be defined before they are invoked, so this subroutine needs to be initialized before we parse arguments below.

printf "\v\t This function runs connectivity of an ROI to every other voxel in the brain based on standardized connectomes."
printf "%15s\n\v\t -h  | --help"
printf "\t Print this usage text"
printf "%15s\n\t -r | --roi <roi.nii.gz, roi.nii or roilist.txt>"
printf "\t Define the regions that will be used as ROIs. To run multiple regions, supply a textfile"
printf "%15s\n\t -f | --fMRI"
printf "\t Calculate connectivity using fMRI. If neither -f nor -d are supplied, connectivity will be done for both fMRI and dMRI."
printf "%15s\n\t -d | --dMRI"
printf "\t Calculate connectivity using dMRI. If neither -f nor -d are supplied, connectivity will be done for both fMRI and dMRI."
printf "%15s\n\t -c | --command <seed/pseed/matrix/pmatrix/seedvox_noram/seedvox_ram>"
printf "%15s\n\t Define what to do with the ROI(s). Defaults to seed which will run connectivity from that seed to the rest of the brain."
printf "%15s\n\t pseed runs connectivity from first seed but partialing out all other seeds (.txt file needed in -r)"
printf "%15s\n\t matrix exports a connectivity matrix between all seeds supplied (.txt file needed in -r)"
printf "%15s\n\t pmatrix exports a partial connectivity matrix between all seeds supplied (.txt file needed in -r)"
printf "%15s\n\t seedvox_noram exports weighted averaged connectivity seeding from each nonzero voxel of the seed. This takes a long time. Use extended queue."
printf "%15s\n\t seedvox_ram exports weighted averaged connectivity seeding from each nonzero voxel of the seed but inverts the whole brain matrix in RAM. Supply max RAM available and use small ROIs."
printf "%15s\n\t -s | --writesinglefiles"
printf "\t (Only) when supplied, single subject output is written, too."
printf "%15s\n\t -o | --outputfolder </path/to/folder/>"
printf "\t Define where to output files to."
printf "%15s\n\t -m | --mask <maskfile.nii/.nii.gz>"
printf "\t Define the area of which to return output to. If left blank, the whole brain will be calculated."

printf "\v\t EXAMPLE : ./connectome.sh -r /autofs/cluster/nimlab/rois/testroi.nii.gz -o /autofs/cluster/nimlab/output/ -f -d -c seed"
printf "\n\t will run fMRI/dMRI based connectivity seeding from the testroi.nii.gz seed and output will be placed in the output directory \n\n"

}

# define defaults here
writesingle="0"
command="seed"
outputfolder=$HOME
maskname="."
doboth=1
dodMRI="0"
dofMRI="0"

[[ "$1" == "" ]] && usage && exit
while [ "$1" != "" ]; do
    case $1 in
        -h  | --help)           shift
                                usage #RUN USAGE HERE
                                exit
                                ;;
        -r | --roi )            shift
                                filename=$1
                                ;;
        -m | --mask )           shift
                                maskname=$1
                                ;;
        -c | --command )        shift
			                    command=$1
                                ;;
        -s | --writesinglefiles )    writesingle="1"
                                ;;
        -o | --outputfolder )   shift
                                outputfolder=$1
                                ;;
        -f | --fMRI )           dofMRI="1" 
                                doboth=0
                                ;;
        -d | --dMRI )           dodMRI="1" 
                                doboth=0
    esac
    shift
done

if [ -z $filename ]
    then

        echo "Please supply a seed filename (.nii/.nii.gz/.txt)!"
        exit
    fi

        if [ $outputfolder == $HOME ]
            then
            # default outputfolder to path from seed:
            outputfolder=$(dirname "${filename}")
        fi

    if [ $doboth == 1 ]
        then
            dofMRI="1"
            dodMRI="1"
    fi

        if [ $dofMRI == "1" ] && [ $dodMRI == "1" ]; then
            doboth=1
        fi
        # setup command here:
        cmd="/home/agh14/lead_dbs/connectomics/mapper/run_cs_conseed.sh /usr/pubsw/common/matlab/8.6 $dofMRI $dodMRI /n/scratch2/agh14/connectomes/ $filename $command $writesingle $outputfolder $maskname"


        filebasename=$(basename "$filename")
        extension="${filebasename##*.}"
        #echo $extension
        if [ $extension == "txt" ] && [ $command == "seed" ]; then # multiple seeds, read in and supply separately.
            while IFS='' read -r line || [[ -n "$line" ]]; do
                if [ $doboth == 1 ]
                    then # split jobs for fMRI and dMRI
                    cmd="/home/agh14/lead_dbs/connectomics/mapper/run_cs_conseed.sh /usr/pubsw/common/matlab/8.6 1 0 /n/scratch2/agh14/connectomes/ $line $command $writesingle $outputfolder $maskname"
                    echo $cmd
                    pbsubmit -q highio -l vmem=30gb -c "$cmd"
                    cmd="/home/agh14/lead_dbs/connectomics/mapper/run_cs_conseed.sh /usr/pubsw/common/matlab/8.6 0 1 /n/scratch2/agh14/connectomes/ $line $command $writesingle $outputfolder $maskname"
                    echo $cmd
                    pbsubmit -q highio -l vmem=30gb -c "$cmd"
                else
                cmd="/home/agh14/lead_dbs/connectomics/mapper/run_cs_conseed.sh /usr/pubsw/common/matlab/8.6 $dofMRI $dodMRI /n/scratch2/agh14/connectomes/ $line $command $writesingle $outputfolder $maskname"
                echo $cmd
                pbsubmit -q highio -l vmem=30gb highio -c "$cmd"
                fi
            done < "$filename"
        else
            if [ $doboth == 1 ]
                then # split jobs for fMRI and dMRI
                    cmd="/home/agh14/lead_dbs/connectomics/mapper/run_cs_conseed.sh /usr/pubsw/common/matlab/8.6 1 0 /n/scratch2/agh14/connectomes/ $filename $command $writesingle $outputfolder $maskname"
                    echo $cmd
                    pbsubmit -q highio -l vmem=30gb -c "$cmd"
                    cmd="/home/agh14/lead_dbs/connectomics/mapper/run_cs_conseed.sh /usr/pubsw/common/matlab/8.6 0 1 /n/scratch2/agh14/connectomes/ $filename $command $writesingle $outputfolder $maskname"
                    echo $cmd
                    pbsubmit -q highio -l vmem=30gb -c "$cmd"
            else
                echo $cmd
                pbsubmit -q highio -l vmem=30gb -c "$cmd"
            fi
        fi






function usage
{
    echo "usage: system_page [[[-r roi ] [-c command]] | [-h help]]"
}
