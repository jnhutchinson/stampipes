###################
# Creates marked BAM file
###################
# REQUIRED VARIABLES
###################
# SAMPLE_NAME=Example_NoIndex_L007
###################
# OPTIONAL VARIABLES
# Note: SAMPLE_NAME is not necessary if BAMFILE and DUP_OUT are specified.
###################
# BAMFILE : The input BAM files (default $(SAMPLE_NAME).sorted.bam)
# DUP_OUT : The output file (default $(SAMPLE_NAME).sorted.bam.picard.dup)
# OUTBAM : The output bam (default: /dev/null)
###################
# REQUIRED MODULES
###################
# module load jdk
# module load picard
###################

BAMFILE ?= $(SAMPLE_NAME).uniques.sorted.bam
DUP_OUT ?= $(SAMPLE_NAME).MarkDuplicates.picard
OUTBAM ?= $(SAMPLE_NAME).uniques.sorted.marked.bam
TMPDIR ?= $(shell pwd)

all : $(OUTBAM)

# Add mate cigar information
$(TMPDIR)/$(SAMPLE_NAME).cigar.bam : $(BAMFILE)
	time picard RevertOriginalBaseQualitiesAndAddMateCigar \
		INPUT=$(BAMFILE) OUTPUT=$(TMPDIR)/$(SAMPLE_NAME).cigar.bam \
		VALIDATION_STRINGENCY=SILENT RESTORE_ORIGINAL_QUALITIES=false SORT_ORDER=coordinate MAX_RECORDS_TO_EXAMINE=0

# Calculate the duplication score of the random sample
$(OUTBAM) : $(TMPDIR)/$(SAMPLE_NAME).cigar.bam
	time picard MarkDuplicatesWithMateCigar INPUT=$(TMPDIR)/$(SAMPLE_NAME).cigar.bam OUTPUT=$(OUTBAM) \
		METRICS_FILE=$(DUP_OUT) ASSUME_SORTED=true MINIMUM_DISTANCE=150 VALIDATION_STRINGENCY=SILENT \
		READ_NAME_REGEX='[a-zA-Z0-9]+:[0-9]+:[a-zA-Z0-9]+:[0-9]+:([0-9]+):([0-9]+):([0-9]+).*'
