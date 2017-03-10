###################
# Creates a duplication score from a sorted BAM file.
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
OUTBAM ?= /dev/null
OUTBAMIND ?= $(OUTBAM).bai
TMPDIR ?= $(shell pwd)

all : $(DUP_OUT) $(OUTBAM) $(OUTBAMIND)

# Get rid of multi-mapping, low quality, and chrM alignments
$(TMPDIR)/$(SAMPLE_NAME).nochrM.bam : $(BAMFILE)
	samtools view -F 512 -b $(BAMFILE) chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr1 chr20 chr21 chr22 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chrX chrY > \
$(TMPDIR)/$(SAMPLE_NAME).nochrM.bam

# Add mate cigar information
$(TMPDIR)/$(SAMPLE_NAME).nochrM.cigar.bam : $(TMPDIR)/$(SAMPLE_NAME).nochrM.bam
	time picard RevertOriginalBaseQualitiesAndAddMateCigar \
		INPUT=$(TMPDIR)/$(SAMPLE_NAME).nochrM.bam OUTPUT=$(TMPDIR)/$(SAMPLE_NAME).nochrM.cigar.bam \
		VALIDATION_STRINGENCY=SILENT RESTORE_ORIGINAL_QUALITIES=false SORT_ORDER=coordinate

# Calculate the duplication score of the random sample
$(DUP_OUT) : $(TMPDIR)/$(SAMPLE_NAME).nochrM.cigar.bam
	time picard UmiAwareMarkDuplicatesWithMateCigar INPUT=$(TMPDIR)/$(SAMPLE_NAME).nochrM.cigar.bam OUTPUT=$(TMPDIR)/$(SAMPLE_NAME).nochrM.cigar.dupe.bam \
		METRICS_FILE=$(DUP_OUT) UMI_TAG_NAME=XD ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT \
		READ_NAME_REGEX='[a-zA-Z0-9]+:[0-9]+:[a-zA-Z0-9]+:[0-9]+:([0-9]+):([0-9]+):([0-9]+).*'

# Get rid of optical/PCR dupes
$(OUTBAM) : $(TMPDIR)/$(SAMPLE_NAME).nochrM.cigar.dupe.bam
	samtools view -b -F 1024 $(TMPDIR)/$(SAMPLE_NAME).nochrM.cigar.dupe.bam > $(OUTBAM)

# Index
$(OUTBAMIND) : $(OUTBAM)
	samtools index $(OUTBAM)