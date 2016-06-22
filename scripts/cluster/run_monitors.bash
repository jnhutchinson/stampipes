#!/bin/bash -l

source ~/.bashrc
source $PYTHON3_ACTIVATE

date +"%Y-%m-%d_%H-%M-%S"

python3 $STAMPIPES/scripts/cluster/monitor_alignments.py
