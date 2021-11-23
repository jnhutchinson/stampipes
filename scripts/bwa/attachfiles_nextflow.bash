# Uploads important files as attachments to objects on the LIMS
# Required environment names:
#  * STAMPIPES
#  * ALIGN_DIR
#  * ALIGNMENT_ID

UPLOAD_SCRIPT=$STAMPIPES/scripts/lims/upload_data.py

function attach_alignment () {
    python3 "$UPLOAD_SCRIPT" \
        --attach_file_contenttype SequencingData.flowcelllanealignment \
        --attach_file_objectid "$ALIGNMENT_ID" \
        "$@"
}
function attach_aln_file () {
    local filename=$1
    local purpose=$2
    local filetype=$3
    attach_alignment \
        --attach_file "$filename"  \
        --attach_file_purpose "$purpose" \
        --attach_file_type "$filetype"
}

attach_alignment --attach_directory "${ALIGN_DIR}" --attach_file_purpose alignment-directory

attach_aln_file marked.cram         all-alignments-bam           cram
attach_aln_file density.bw          density-bigwig-windowed      bigwig
attach_aln_file density.bed.starch  density-bed-starch-windowed  starch
