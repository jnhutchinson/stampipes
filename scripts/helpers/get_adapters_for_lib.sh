#!/bin/bash

source "$(dirname "$0")"/lims/api_functions.sh

echo -e "#Lane_ID\tLibNum\tBarcode1\tBarcode2\tP7\tP5\tlane_dir\tmsg"

for lib_num in "$@" ; do

  lims_get_all "flowcell_lane/?library__number=$lib_num" \
    | jq .id \
    | while read -r lane_id ; do
      echo -ne "$lane_id\tLN$lib_num\t"
      lims_get "flowcell_lane/$lane_id/processing_information/" \
        | jq -r -c '.libraries[0] | [
          .barcode1.label_id,
          if .barcode2 then
            .barcode2.label_id
          else
            "."
          end,

          .barcode1.adapter7,
          if .barcode2 then
            .barcode2.adapter5_reverse_complement
          else 
            .barcode1.adapter5
          end,
          .directory
          ] | join ("\t")' \
      #| jq '{"barcode1": .libraries[0].barcode1, "barcode2": .libraries[0].barcode2}'
  done

done | while read -r lane lib bc1 bc2 p7 p5 dir ; do
  echo -n $lane $lib $bc1 $bc2 $p7 $p5 ""
  
  adapter_file=$(ls "$dir"/align*/*adapters.txt 2>/dev/null | sort | head -n1 )
  if [[ ! -s "$adapter_file" ]] ; then
    echo "no_adapter_file,"
    continue
  fi

  disk_p7=$(awk '$1=="P7" {print $2}' "$adapter_file")
  disk_p5=$(awk '$1=="P5" {print $2}' "$adapter_file")

  err=
  if [[ "$p7" != "$disk_p7" ]] ; then
    err=p7_differ,
  fi
  if [[ "$p5" != "$disk_p5" ]] ; then
    err=${err}p5_differ,
  fi
  if [[ -n "$err" ]] ; then
    echo "$err"
  else
    echo ":)"
  fi

done \
  | sed 's/\s\+/	/g'
#https://lims.altiusinstitute.org/api/flowcell_lane/58940/processing_information/
