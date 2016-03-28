#!/bin/bash

die(){
  echo $@ >&2
  exit 1
}
if [ $# -eq 0 ] ; then
  die "Usage: $0 lib_num [lib_nums...]"
fi


if [[ ! -d $CUSTOMER_DISTRIBUTION ]] ;then
  die "\$CUSTOMER_DISTRIBUTION envvar must be set and exist"
fi

source "$STAMPIPES/scripts/lims/api_functions.sh"

fastq_purpose=$(lims_get_all "file_purpose/?slug=fastqc-directory" | jq -r .id)

for libnum in "$@"; do

  libnum=${libnum/LN/}

  lane_ids=$(lims_get_all "flowcell_lane/?library__number=$libnum&warn=False&failed=False&flowcell__failed=False" | jq -r .id )

  echo "# LN$libnum - $(wc -w <<< $lane_ids) lanes"

  for lane_id in $lane_ids; do
    echo " # lane $lane_id"
    lane_info=$(lims_get "flowcell_lane/$lane_id/")
    content_type=$(jq -r .object_content_type <<< $lane_info)
    lane_path=$(lims_get_all "directory/?purpose=$fastq_purpose&content_type=$content_type&object_id=$lane_id" | jq -r .path )
    flowcell="FC$(jq -r .flowcell_label <<< $lane_info )"
    project_url=$(jq -r .project <<< $lane_info )
    cust_dir="$(lims_get_by_url $project_url | jq -r .distribution_directory)"
    dist_dir="$CUSTOMER_DISTRIBUTION/$cust_dir"

    if [[ -z "$cust_dir" || ! -d $dist_dir ]] ;then
      die "Error: Distribution directory '$dist_dir' for project '$project_url' - library '$libnum,' lane '$lane_id', does not exist" >&2
    fi
    if [[ ! -d $lane_path ]] ;then
      die "Error: Lane path '$lane_path' doesn't exist for library 'LN$libnum', lane '$lane_id'"
    fi

    echo "   " mkdir -p "$dist_dir/$flowcell"
    echo "   " ln -s "$lane_path" "$dist_dir/$flowcell"
  done
  echo

done
