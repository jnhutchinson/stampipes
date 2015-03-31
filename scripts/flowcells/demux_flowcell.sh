#!/bin/bash
set -e

usage(){
cat << EOF
usage: $0 [-f flowcell_label]

Setup analysis for a flowcell

OPTIONS:
-h        Show this message
-i        Unaligned input directory
-o        Unaligned output directory
-p        processing.json
EOF
}

check_folder(){
  if [ ! -d "$1" ] ; then
    echo "Directory $1 does not exist, aborting." >&2
    exit 1
  fi
}

demux_script="$STAMPIPES/scripts/flowcells/demux_fastq.py"
indir=
outdir=
processing=
mismatches=0
while getopts ":hi:o:p:m:" opt ; do
  case $opt in
    h)
      usage
      exit 0
      ;;
    i)
      indir=$(readlink -f "$OPTARG")
      ;;
    o)
      outdir="$OPTARG"
      ;;
    p)
      processing=$(readlink -f "$OPTARG")
      ;;
    m)
      mismatches="$OPTARG"
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

if [ -z "$processing" ] || [ -z "$indir" ] || [ -z "$outdir" ] ; then
  usage >&2
  exit 1
fi

if [ ! -s "$processing" ] ; then
  echo "$processing is not a non-empty file" >&2
  exit 1
fi

inputfiles=($(find "$indir" -name '*Undetermined_*fastq.gz'))

for i in "${inputfiles[@]}" ; do
  lane=$( sed 's/.*_L\(00.\)_.*/\1/' <(basename "$i" ))

  qsub -cwd -V -q all.q -N .dmx$(basename $i) <<__DEMUX__
    python "$demux_script"        \
      --autosuffix                \
      --processing "$processing"  \
      --outdir "$outdir"          \
      --mismatches "$mismatches"  \
      --lane "$lane"              \
      "$i"
__DEMUX__

done
