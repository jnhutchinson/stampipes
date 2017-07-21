###################
# These variables must be passed in or set for the makefile to work.
###################
# FASTQ_FILE=/path/to/file.fastq.gz
# FASTQC_FILE=/path/to/file_fastqc.zip
###################
###################

# Use the default PATH version of each of these programs unless overridden
FASTQC ?= fastqc

# Use python3
PYTHON ?= python3

# Use eight threads by default; each thread uses 250MB of memory, so we are staying
# well within qsub memory limits
FASTQC_OPTIONS ?= -t 8 --noextract --nogroup --casava --limits $(STAMPIPES)/config/fastqc/limits.txt

ORIGINAL_FASTQC_OUTFILE ?= $(basename $(basename $(FASTQ_FILE) ))_fastqc.zip

FASTQC_FILE ?= $(ORIGINAL_FASTQC_OUTFILE)

OUTDIR ?= $(shell pwd)

STAMPIPES ?= ~/stampipes

.PHONY : info

all : info fastqc

fastqc : $(FASTQC_FILE)
	
info :
	@echo "------"
	@echo "PROGRAM VERSIONS"
	@echo "------"
	@echo "FastQC:"
	@$(FASTQC) -v
	@echo "------"
	@echo "METADATA"
	@echo "------"
	@echo "FILE: " $(FASTQ_FILE)
	@echo "------"
	@echo "PROGRAM OPTIONS"
	@echo "------"
	@echo "FASTQC_OPTIONS: " $(FASTQC_OPTIONS)
	@echo "------"

$(FASTQC_FILE) : $(FASTQ_FILE)
	@echo FastQ file: $(FASTQ_FILE)
	time $(FASTQC) -o $(OUTDIR) $(FASTQC_OPTIONS) $(FASTQ_FILE) && echo FastQC stats >&2
ifneq ($(FASTQC_FILE),$(ORIGINAL_FASTQC_OUTFILE))
	@echo Moving to $(FASTQC_FILE) from $(ORIGINAL_FASTQC_OUTFILE)
	mv $(ORIGINAL_FASTQC_OUTFILE) $(FASTQC_FILE)
endif
