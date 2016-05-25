# !/bin/bash

# create project and genome folders with requisite files if they do not exist

PROJECT=$1
GENOME=$2
OUTPUTDIR=/net/seq/data/flowcells/trackhubs/

# if project doesn't exist, create it and its necessary files
if [ ! -d $OUTPUTDIR$PROJECT ];
then
    mkdir $OUTPUTDIR$PROJECT
    echo "hub master_hub_$PROJECT" > $OUTPUTDIR$PROJECT/hub.txt
    echo "shortLabel MasterHub_$PROJECT" >> $OUTPUTDIR$PROJECT/hub.txt
    echo "longLabel list of all hubs..." >> $OUTPUTDIR$PROJECT/hub.txt
    echo "genomesFile genomes.txt" >> $OUTPUTDIR$PROJECT/hub.txt
    echo "email placeholder@altiusinstitute.org" >> $OUTPUTDIR$PROJECT/hub.txt

    echo "" > $OUTPUTDIR$PROJECT/genomes.txt
fi

# if genome for the project doesn't exist, create it, create its necessary folders, and add it to the genome list
if [ ! -d $OUTPUTDIR$PROJECT/$GENOME ];
then
    mkdir $OUTPUTDIR$PROJECT/$GENOME
    mkdir $OUTPUTDIR$PROJECT/$GENOME/backup_master_tracks
    mkdir $OUTPUTDIR$PROJECT/$GENOME/tracks_by_flowcell
    echo "genome $GENOME" >> $OUTPUTDIR$PROJECT/genomes.txt
    echo "trackDb $GENOME/master_track.txt" >> $OUTPUTDIR$PROJECT/genomes.txt
    echo "" >> $OUTPUTDIR$PROJECT/genomes.txt
fi