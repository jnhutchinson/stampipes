###################
# These variables must be passed in or set for the makefile to work.  If the genome's
# FAI file is not at $(BWAINDEX).fai, then it must also be specified under FAI.
###################
# SAMPLE_NAME := Sample_NoIndex_L007
# BWAINDEX := /path/to/genome
###################
# REQUIRED MODULES
###################
# module load java
# module load samtools
# module load python
###################

SHELL = bash

# Use the default PATH version of each of these programs unless overridden
SAMTOOLS ?= samtools

FAI ?= $(BWAINDEX).fai

READ_LENGTH ?= 36

MIN_MAP_QUALITY ?= 30

SAMTOOL_OPTIONS ?= -F 0x12 -q $(MIN_MAP_QUALITY)

# Ideally this will be set to something else set in the environment or on the command line
TMPDIR ?= .

# where our results files go
OUTDIR ?= .

all : info metrics uniques

info : 
	@echo "------"
	@echo "PROGRAM VERSIONS"
	@echo "------"
	@echo "samtools:"
	@$(SAMTOOLS) 2>&1 | grep "Version"
	@echo "------"
	@echo "METADATA"
	@echo "------"
	@echo "SAMPLE_NAME: " $(SAMPLE_NAME)
	@echo "BWAINDEX: " $(BWAINDEX)
	@echo "FAI: " $(FAI)
	@echo "READ LENGTH: " $(READ_LENGTH)

metrics : $(OUTDIR)/$(SAMPLE_NAME).CollectInsertSizeMetrics.picard

uniques : $(OUTDIR)/$(SAMPLE_NAME).uniques.sorted.bam.bai

# Sometimes this will report errors about a read not mapping that should have a mapq of 0
# See this for more info: http://seqanswers.com/forums/showthread.php?t=4246
$(OUTDIR)/$(SAMPLE_NAME).CollectInsertSizeMetrics.picard : $(OUTDIR)/$(SAMPLE_NAME).sorted.bam
	time java -Xmx1000m -jar `which CollectInsertSizeMetrics.jar` INPUT=$^ OUTPUT=$@ \
                HISTOGRAM_FILE=$(SAMPLE_NAME).CollectInsertSizeMetrics.pdf \
                VALIDATION_STRINGENCY=LENIENT \
                ASSUME_SORTED=true && echo Picard stats >&2

# Index uniquely mapping reads 
$(OUTDIR)/$(SAMPLE_NAME).uniques.sorted.bam.bai : $(OUTDIR)/$(SAMPLE_NAME).uniques.sorted.bam
	time $(SAMTOOLS) index $^
                
# Sorted uniquely mapping reads BAM
$(OUTDIR)/$(SAMPLE_NAME).uniques.sorted.bam : $(OUTDIR)/$(SAMPLE_NAME).sorted.bam
	time $(SAMTOOLS) view $(SAMTOOL_OPTIONS) $^ | $(SAMTOOLS) view -bSt $(FAI) - > $@

# Index sorted BAM file
$(OUTDIR)/$(SAMPLE_NAME).sorted.bam.bai : $(OUTDIR)/$(SAMPLE_NAME).sorted.bam
	time $(SAMTOOLS) index $^	
