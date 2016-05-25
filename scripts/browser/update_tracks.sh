# !/bin/bash

# backup and update master tracks for given project/genome

PROJECT=$1
GENOME=$2
OUTPUTDIR=/net/seq/data/flowcells/trackhubs/

timestamp=$(date +%s)
# list of track names to create master list from, the 'master_track_list.txt' is updated in this script to contain all file names
# but a seperate list of desired names can be supplied
buildfile="master_track_list.txt"

if [ ! -d $OUTPUTDIR$PROJECT/$GENOME ];
then
    exit;
fi

# backup master track
if [ -f $OUTPUTDIR$PROJECT/$GENOME/master_track.txt ];
then
    mv $OUTPUTDIR$PROJECT/$GENOME/master_track.txt $OUTPUTDIR$PROJECT/$GENOME/backup_master_tracks/backup_master_track_$timestamp\.txt
fi

# create new build list of all names
> $OUTPUTDIR$PROJECT/$GENOME/master_track_list.txt
cd $OUTPUTDIR$PROJECT/$GENOME/tracks_by_flowcell
for filename in *
do
    echo "$filename" >> $OUTPUTDIR$PROJECT/$GENOME/master_track_list.txt
done

# create new master track given the build list
> $OUTPUTDIR$PROJECT/$GENOME/master_track.txt
while IFS= read line
do
    cat $OUTPUTDIR$PROJECT/$GENOME/tracks_by_flowcell/$line >> $OUTPUTDIR$PROJECT/$GENOME/master_track.txt
done < $OUTPUTDIR$PROJECT/$GENOME/$buildfile