#!/bin/bash

source "$(dirname "$0")/api_functions.sh"
dryrun=
while getopts 'n' opt; do
  case "$opt" in
    n)
      dryrun='true'
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 2
  esac
done
shift $((OPTIND-1))

lane_id="$1"
r1_file="${2:-}"
r2_file="${3:-}"

root="/net/seq/data/flowcells"
#root=/net/monarch/vol2/tag/stamlab/flowcells

authoritative_dir="/net/seq/data/data-release/encode/DNase/ /net/seq/data/data-release/encode/DGF/"

if [ -z "$lane_id" ] ; then
  echo "Usage: $0 lane_id"
  exit 2
fi

# Tries to upload, returns if no file
# Exits if successful
function try_upload(){
  if [[ -z "$r1_file" ]] ; then
    return
  fi
  cmd=python3
  if [[ -n $dryrun ]] ; then
    cmd=echo
  fi
  #python3 "$STAMPIPES/scripts/lims/upload_data.py" --attach_file_contenttype SequencingData.flowcelllane --attach_file_objectid "$lane_id" --attach_directory "$dir" --attach_file_purpose fastq-directory
  if [ -s "$r1_file" ] ; then
    "$cmd" "$STAMPIPES/scripts/lims/upload_data.py" --attach_file_contenttype SequencingData.flowcelllane --attach_file_objectid "$lane_id" --attach_file_type gzipped-fastq --attach_file_purpose r1-fastq --attach_file "$r1_file" #--skip_md5_check
  fi
  if [ -s "$r2_file" ] ; then
    "$cmd" "$STAMPIPES/scripts/lims/upload_data.py" --attach_file_contenttype SequencingData.flowcelllane --attach_file_objectid "$lane_id" --attach_file_type gzipped-fastq --attach_file_purpose r2-fastq --attach_file "$r2_file" #--skip_md5_check
  fi
  exit
}

# Check to see if files exist already
echo "# Lane: $lane_id"
info=$(lims_get flowcell_lane/$lane_id/)
bin=$(jq -r .bin <<< "$info")
lane=$(jq -r .lane <<< "$info")
flowcell=$(jq -r .flowcell_label <<< "$info")
project=$(jq -r .project_label <<< "$info")
sample=$(jq -r .samplesheet_name <<< "$info")

file_count=$(lims_get "file/?object_id=$lane_id&content_type=40" | jq -r .count)
if [ "$file_count" -gt 0 ] ; then
  echo "# Already has files"
  exit 0
fi

info2=$(lims_get flowcell_lane/$lane_id/processing_information/)
fcdir=$(jq -r .flowcell.directory <<< "$info2")


# Check canonical flowcell directory
if [  "$project" == "null" ] ; then
  project=*
fi

# Find flowcell directory
if [ "$fcdir" == "null" ] ; then
  fcdir=$(ls -d $root/FC$flowcell*)
  if [ -z "$fcdir" ] ; then
    exit 1
  fi
fi
echo "# dir(s): $fcdir"


# Find sample directory
dir=$(ls -d $fcdir/Project_$project/Sample_$sample* 2>/dev/null)
if [ ! -d "$dir" ] ; then
  echo "Trying w/o sublibrary"
  sample=$(sed 's/[A-Z]$//' <<< $sample)
  dir=$(ls -d $fcdir/Project_$project/Sample_$sample* 2>/dev/null)
fi
echo "Dir is $dir"

# Look for collated files
if [ -d "$dir" ] && [ -n "$(ls "$dir"/*gz 2>/dev/null)" ] ; then
  r1_file=$(ls "$dir/"*"$sample"*_L00"$lane"*_R1*.fastq.gz)
  r2_file=$(ls "$dir/"*"$sample"*_L00"$lane"*_R2*.fastq.gz)

  if [ "$(wc -w <<< $r1_file)" -gt 1 ] ; then
    echo "! Need collation: $lane_id $dir"
    exit
  fi
fi

try_upload

# Look for bin file
if [ "$bin" != "null" ] ; then
  binstub=$(printf "%03d" "$bin")
  dir=$(ls -d $fcdir/$binstub)
  r1_file="$( ls $root/FC$flowcell*/$binstub/s_${lane}_sequence.txt.gz | tail -n1)"
  try_upload
fi

# Look for sequence files
r1_file="$(ls "$fcdir/s_${lane}_1_sequence.txt.gz" -S 2>/dev/null | head -n1)"
r2_file="$(ls "$fcdir/s_${lane}_2_sequence.txt.gz" -S 2>/dev/null | head -n1)"
echo "r1: $r1_file"
if [ -z "$r1_file" ] ; then
  r1_file=$(ls $fcdir/s_${lane}_sequence.txt.gz -S | head -n1)
fi

try_upload

# Look in archival directories
echo '#' $lane_id - $sample - $flowcell - $lane
pattern="*$sample*${flowcell}_$lane*fastq.gz*"
actual_files=$(find $authoritative_dir -name "$pattern" | grep -v ".bak" | grep -v "old")
#actual_files=$(find $authoritative_dir -name "$pattern" )
echo "# Looking for $pattern"
if [ -n "$actual_files" ] ; then
 for f in $actual_files ; do
   echo $f  ":" $(readlink -f $f)
 done
 if [[ $(wc -l <<< "$actual_files") -eq 1 ]] ; then
   r1_file=$(readlink -f "$actual_files")
   echo "R1! $r1_file"
 else
   echo "Too many!"
 fi
fi

try_upload


echo "FAILURE: Did not find a fastq file for $lane_id :("
