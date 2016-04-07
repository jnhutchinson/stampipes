This is a focused collection of makefiles and scripts to for analysis on Illumina sample data. While this probably can't be directly used by someone not in our lab, it might be a useful reference for others.

Set up
-------------

In order to set up these scripts, checkout the repository.  There are some environment variables to set:

* STAMPIPES -- full path to the installation directory
* LIMS_API_URL -- URL to the LIMS API
* LIMS_API_TOKEN -- Auth to the LIMS API
* MODULELOAD -- Full path to the GNU modules init script
* HOTSPOT_DIR -- Path to HOTSPOT, defaults to `~/hotspot/hotspot-distr` - download at [https://github.com/rthurman/hotspot]
* PYTHON3_ACTIVATE -- A virtualenv activation script for python3.

A python version that includes pip needs to be installed for proper python package version recording.

Modules
-------------

The current modules to load for working the pipeline are:

    module load bedops/2.4.15
    module load bedtools/2.16.2
    module load bowtie/1.0.0
    module load bwa/0.7.12
    module load coreutils/8.9
    module load cufflinks/2.2.1
    module load FastQC/0.11.3
    module load gcc/4.7.2
    module load git/2.3.3
    module load java/jdk1.7.0_05
    module load picard/1.118
    module load R/3.1.0
    module load samtools/1.2
    module load tophat/2.0.13

If $MODULELOAD is set, then the processing scripts should load the required modules automatically.
