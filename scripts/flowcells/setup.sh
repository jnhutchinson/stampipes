#!/bin/bash

# Dependencies
source $MODULELOAD
module load python/2.7.3

#########
# Options
#########

NODES=40
THREADS=8

usage(){
cat << EOF
usage: $0 [-f flowcell_label]

Setup analysis for a flowcell

OPTIONS:
-h        Show this message
-v        Verbose
-f        Flowcell Label
EOF
}

while getopts ":hvf:" opt ; do
    case $opt in
    h)
        usage
        exit 0
        ;;
    v)
        verbose=true
        ;;
    f)
        flowcell="$OPTARG"
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

if [ -z "$flowcell" ] ; then
    echo "No flowcell label specified"
    flowcell=$(basename $(pwd) | cut -f4 -d_ | cut -c2-6)
    echo "Guessing $flowcell..."
    #usage >&2
    #exit 1
fi

#######################
# Samplesheet functions
#######################

make_hiseq_samplesheet(){
echo "FCID,Lane,SampleID,SampleRef,Index,Description,Control,Recipe,Operator,SampleProject"

# Escape double quotes, as we need to interpolate $flowcell
# ( X | tostring) syntax is for non-string fields
# Default values (if field is false or null) come after //
jq -r ".libraries[] | [
 \"FC$flowcell\",
 (.lane | tostring),
 .samplesheet_name,
 .alignments[0].genome_index // \"contam\",
 .barcode_index,
 .cell_type                  // \"None\"  ,
 \"N\",
 .assay                      // \"N/A\"   ,
 \"orders\",
 .project
 ] | join(\",\") " "$json"
}

make_nextseq_samplesheet(){
    name=Stamlab
    date=$(date '+%m/%d/%Y')
cat <<__SHEET__
[Header]
Investigator Name,$name
Project Name,$name
Experiment Name,$name
Date,$date
Workflow,GenerateFASTQ

[Settings]

[Data]
SampleID,SampleName,index,index2
__SHEET__

# This bit of cryptic magic generates the samplesheet part.
jq -r '.libraries[] | [.samplesheet_name,.samplesheet_name,.barcode_index,""] | join(",") ' "$json"

}

########
# Main #
########

json="processing.json"
illumina_dir=$(pwd)

link_command="#no linking to do"

# Get and read the processing script
python "$STAMPIPES/scripts/lims/get_processing.py" -f "$flowcell" -o "$json"
run_type=$(     jq -r '.flowcell.run_type'          "$json" )
analysis_dir=$( jq -r '.alignment_group.directory'  "$json" )
mask=$(         jq -r '.alignment_group.bases_mask' "$json" )
run_type=$(     jq -r '.flowcell.run_type'          "$json" )
barcodes=$(     jq -r '.libraries[].barcode_index'  "$json" )

mismatches=$( "$STAMPIPES/scripts/max_mismatch.py" $barcodes )

case $run_type in
"NextSeq 500")
    echo "Regular NextSeq 500 run detected"
    parallel_env="-pe threads $THREADS"
    link_command="python $STAMPIPES/scripts/flowcells/link_nextseq.py -i fastq -o ."
    samplesheet="SampleSheet.csv"
    fastq_dir="$illumina_dir/fastq"  # Lack of trailing slash is important for rsync!
    make_nextseq_samplesheet > SampleSheet.csv

    # The quadruple-backslash syntax on this is messy and gross.
    # It works, though, and the output is readable.
    read -d '' unaligned_command  << _U_
    bcl2fastq \\\\
      --input-dir "${illumina_dir}/Data/Intensities/BaseCalls" \\\\
      --use-bases-mask "$mask" \\\\
      --output-dir "$fastq_dir" \\\\
      --with-failed-reads \\\\
      --barcode-mismatches "$mismatches"
_U_
    ;;
    #TODO: Add HISEQ V3 on hiseq 2500 (rapid run mode)
"HISEQ V4")
    echo "Regular HiSeq 2500 run detected"
    parallel_env=""
    link_command='#no linking to do'
    samplesheet=$(pwd)/Data/Intensities/BaseCalls/SampleSheet.csv
    mkdir -p $(basename "$samplesheet")
    make_hiseq_samplesheet > "$samplesheet"
    fastq_dir="$illumina_dir/Unaligned/"  # Trailing slash is important for rsync!

    read -d '' unaligned_command <<_U_
    if [ ! -e "$fastq_dir" ] ; then
            configureBclToFastq.pl \\\\
              --mismatches "$mismatches" \\\\
              --output-dir "$fastq_dir" \\\\
              --fastq-cluster-count 16000000 \\\\
              --with-failed-reads --sample-sheet $samplesheet \\\\
              --use-bases-mask "$mask"  \\\\
              --input-dir "$illumina_dir/Data/Intensities/BaseCalls"
    fi

    cd "$fastq_dir"
    qmake -now no -cwd -q all.q -V -- -j "$NODES"
_U_
    ;;
\?)
    echo "Unrecognized sequencer $sequencer"
    exit 1
    ;;
esac


# The final script is below:
cat > run_bcl2fastq.sh <<__BCL2FASTQ__
#!/bin/bash

source $MODULELOAD
module load bcl2fastq/1.8.4
module load bcl2fastq2/2.15.0.4
module load python/2.7.3

while [ ! -e "$illumina_dir/RTAComplete.txt" ] ; do sleep 60 ; done

qsub -cwd -N "u-$flowcell" $parallel_env  -V -S /bin/bash  <<__FASTQ__

set -x -e -o pipefail

cd "$illumina_dir"

$unaligned_command

__FASTQ__

qsub -cwd -N "c-$flowcell" -hold_jid "u-$flowcell" -V -S /bin/bash <<__COPY__
mkdir -p "$analysis_dir"

rsync -avP "$fastq_dir" "$analysis_dir/"
rsync -avP "$illumina_dir/InterOp" "$analysis_dir/"
rsync -avP "$illumina_dir/RunInfo.xml" "$analysis_dir/"
rsync -avP "$samplesheet" "$analysis_dir"

cd "$analysis_dir"

python "$STAMPIPES/scripts/lims/get_processing.py" -f "$flowcell"

$link_command

python "$STAMPIPES/scripts/create_processing.py"

bash run.bash
__COPY__

__BCL2FASTQ__

if [ -e "RTAComplete.txt" ] ; then
    echo -e "Setup complete. To kick everything off, type:\n\nbash run_bcl2fastq.sh"
else
    echo -e "Setup complete, sequencing still in progress. To queue everything up, type:\n\nnohup bash run_bcl2fastq.sh &"
fi
