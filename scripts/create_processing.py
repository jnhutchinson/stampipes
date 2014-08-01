import json
import os

PROCESSING_FILE = 'processing.json'
QSUB_SCRIPT = 'run.bash'
PROCESS_SCRIPT = 'run.bash'

STAMPIPES = os.getenv('STAMPIPES', '~/stampipes')
base_script = open("%s/processes/bwa/process_bwa_paired_trimmed.bash" % STAMPIPES, "r").read()

p = json.loads(open(PROCESSING_FILE, 'r').read())

processing_scripts = dict()

def add_script(script_file, samplesheet_name, priority):
    global processing_scripts

    if not priority in processing_scripts:
        processing_scripts[priority] = list()

    processing_scripts[priority].append((samplesheet_name, script_file))

def run_scripts():
    global processing_scripts

    outfile = open(QSUB_SCRIPT, 'w')

    for priority in sorted(processing_scripts.keys(), reverse=True):
        outfile.write("# Priority %s\n" % str(priority))
        for samplesheet_name, script_file in processing_scripts[priority]:
            outfile.write("cd %s && " % os.path.dirname(script_file))
            outfile.write("qsub -N proc%s -cwd -V -S /bin/bash %s\n\n" % (samplesheet_name, PROCESS_SCRIPT))

    outfile.close()

def create_script(lane):
    global p

    if not lane["alignments"]:
        print "# Skipping %s" % lane['sample']
        return

    alignment = lane["alignments"][0]

    if not alignment['aligner']:
        print "# Skipping %s" % lane['sample']
        return
    print "# Aligning %s" % lane['sample']

    script_file = os.path.join(p['alignment_group']['directory'], "Project_%s" % lane['project'], "Sample_%s" % lane['samplesheet_name'], "run.bash")
    print script_file

    outfile = open(script_file, 'w')
    outfile.write("set -e -o pipefail\n")
    outfile.write("export SAMPLE_NAME=%s\n" % alignment['sample_name'])
    outfile.write("export BWAINDEX=%s\n" % alignment['genome_index_location'])
    outfile.write("export GENOME=%s\n" % alignment['genome_index'])
    outfile.write("export ASSAY=%s\n" % lane['assay'])
    outfile.write("export READLENGTH=%s\n" % p['flowcell']['read_length'])
    if p['flowcell']['paired_end']:
        outfile.write("export PAIRED=True\n")
    outfile.write("export FLOWCELL_LANE_ID=%s\n" % lane['id'])
    outfile.write("export ALIGNMENT_ID=%s\n" % alignment['id'])
    outfile.write("export FLOWCELL=%s\n" % p['flowcell']['label'])
    outfile.write("\n")
    outfile.write(base_script)
    outfile.close()

    add_script(script_file, lane['samplesheet_name'], alignment['priority'])

for lane in p['libraries']:
    create_script(lane)

run_scripts()
