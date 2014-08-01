###################
# These variables must be passed in or set for the makefile to work.
###################
# SAMPLE_NAME=Example_NoIndex_L007
###################
# If the data is paired, then indicate this
###################
# PAIRED=True
###################
# If uploading data to the LIMS, then indicate the LIMS_API_TOKEN to use and FLOWCELL uploading name
# LIMS_API_URL=http://lims.com
# LIMS_API_TOKEN=aethua87ao997i9a9u
# FLOWCELL=CU8TE
###################

# Use the default PATH version of each of these programs unless overridden
FASTQC ?= fastqc

# Use four threads by default; each thread uses 250MB of memory, so we are staying
# well within qsub memory limits
FASTQC_OPTIONS ?= -t 4

INDIR ?= $(shell pwd)

STAMPIPES ?= ~/stampipes

ifdef PAIRED
 TARGETS ?= $(SAMPLE_NAME)_R1_fastqc $(SAMPLE_NAME)_R2_fastqc
else
 TARGETS ?= $(SAMPLE_NAME)_R1_fastqc 
endif 

.PHONY : upload fastqc info

all : info fastqc upload

fastqc : $(TARGETS) 
	@echo "Targets: " $(TARGETS)
	
info :
	@echo "------"
	@echo "PROGRAM VERSIONS"
	@echo "------"
	@echo "FastQC:"
	@$(FASTQC) -v
	@echo "------"
	@echo "METADATA"
	@echo "------"
	@echo "SAMPLE_NAME: " $(SAMPLE_NAME)
	@echo "------"
	@echo "PROGRAM OPTIONS"
	@echo "------"
	@echo "FASTQC_OPTIONS: " $(FASTQC_OPTIONS)
	@echo "------"

# Make fastqc results
# Also move the folder if we only have one file pair so we don't have the extra 001
$(INDIR)/$(SAMPLE_NAME)_%_fastqc :
	time $(FASTQC) $(FASTQC_OPTIONS) $(INDIR)/$(SAMPLE_NAME)_$*_???.fastq.gz --casava && echo FastQC stats >&2
ifeq ($(words $(shell ls $(INDIR)/$(SAMPLE_NAME)_R1_???.fastq.gz | wc -l)), 1)
	mv $(INDIR)/$(SAMPLE_NAME)_$*_001_fastqc $(INDIR)/$(SAMPLE_NAME)_$*_fastqc
	mv $(INDIR)/$(SAMPLE_NAME)_$*_001_fastqc.zip $(INDIR)/$(SAMPLE_NAME)_$*_fastqc.zip
endif

ifdef LIMS_API_TOKEN
upload : $(TARGETS)
	@echo "Uploading FastQC data"
ifdef PAIRED
        python $(STAMPIPES)/scripts/lims/upload_data.py -a $(LIMS_API_URL) -t $(LIMS_API_TOKEN) -f $(FLOWCELL) \
	  --flowcell_lane_id=$(FLOWCELL_LANE_ID) --fastqcfile $(INDIR)/$(SAMPLE_NAME)_R1_fastqc --fastqcfile $(INDIR)/$(SAMPLE_NAME)_R2_fastqc
else
        python $(STAMPIPES)/scripts/lims/upload_data.py -a $(LIMS_API_URL) -t $(LIMS_API_TOKEN) -f $(FLOWCELL) \
	  --flowcell_lane_id=$(FLOWCELL_LANE_ID) --fastqcfile $(INDIR)/$(SAMPLE_NAME)_R1_fastqc
endif
endif
