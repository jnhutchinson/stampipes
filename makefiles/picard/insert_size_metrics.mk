###################
# These variables must be passed in or set for the makefile to work.
###################
# SAMPLE_NAME := Sample_NoIndex_L007
###################
# REQUIRED MODULES
###################
# module load jdk
# module load picard
###################

BAMFILE ?= $(SAMPLE_NAME).uniques.sorted.bam
INSERTMETRICS ?= $(SAMPLE_NAME).CollectInsertSizeMetrics.picard
OUTBAM ?= /dev/null
TMPDIR ?= $(shell pwd)

all: $(INSERTMETRICS)

# Sometimes this will report errors about a read not mapping that should have a mapq of 0
# See this for more info: http://seqanswers.com/forums/showthread.php?t=4246
$(INSERTMETRICS) : $(BAMFILE)
	time picard CollectInsertSizeMetrics INPUT=$(BAMFILE) OUTPUT=$(INSERTMETRICS) \
		HISTOGRAM_FILE=$(INSERTMETRICS).pdf \
		VALIDATION_STRINGENCY=LENIENT \
		ASSUME_SORTED=true && echo Picard stats >&2
