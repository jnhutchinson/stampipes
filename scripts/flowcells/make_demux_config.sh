#!/bin/bash
set -e

usage(){
cat << EOF
usage: $0 [-f flowcell_label]

Setup analysis for a flowcell

OPTIONS:
-h        Show this message
-l        Lane
-u        Flowcell Unaligned directory
-s        Samplesheet.csv
EOF
}

check_folder(){
  if [ ! -d "$1" ] ; then
    echo "Directory $1 does not exist, aborting." >&2
    exit 1
  fi
}

demux_script="$STAMPIPES/scripts/flowcells/demux_fastq.py"
demux_config="demux.config"  # This'll be in $unaligneddir

unaligneddir=
samplesheet=
lane=
while getopts ":hs:u:l:" opt ; do
    case $opt in
    h)
        usage
        exit 0
        ;;
    u)
        unaligneddir="$OPTARG"
        ;;
    s)
        samplesheet="$(readlink -f $OPTARG)"
        ;;
    l)
        lane="$OPTARG"
        ;;
    :)
        echo "Option -$OPTARG requires an argument" >&2
        usage >&2
        exit 1
        ;;
    \?)
        echo "Invalid option: -$OPTARG" >&2
        usage >&2
        exit 1
        ;;
    esac

done


if [ -z "$unaligneddir" ] || [ -z "$samplesheet" ] || [ -z "$lane" ] ; then
	usage >&2
	exit 1
fi

#TODO: Switch between these intelligently

if [ ! -s "$samplesheet" ] ; then
	echo "$samplesheet is not a non-empty file" >&2
	exit 1
fi


if grep -q '^FCID,Lane,SampleID,' "$samplesheet" ; then
	# Hiseq format
  destination="$unaligneddir/Undetermined_indices/Sample_lane$lane"
  check_folder "$destination"
  cd "$destination"
	# TODO: What do we do about dual-index names??
	awk -F, -v "l=$lane" \
		  'BEGIN{OFS="\t"} l == $2 { print $3 "_" $5, $5 }' \
		  "$samplesheet" \
		> "$demux_config"
else
  # NextSeq Format
  destination="$unaligneddir"
  check_folder "$destination"
  cd "$destination"
  awk -F, \
    'BEGIN{OFS="\t";seqline=0}  seqline>0 {print $1 "_S" seqline "_" $3, $3; seqline++} /SampleID,SampleName/ {seqline=1}' \
    "$samplesheet" \
    > "$demux_config"
fi


for i in *Undetermined_*fastq.gz ; do
	echo python "$demux_script" --autosuffix $i --barcodes "$demux_config" \
		| qsub -cwd -V -q all.q -N .dmx$i

done
