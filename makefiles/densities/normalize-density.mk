###################
# Turn a density starch and bigwig into normalized versions.
# 
###################
# SAMPLE_NAME=Example_NoIndex_L007
# BWAINDEX=/path/to/genome/hg19/hg19
# GENOME=hg19
# READLENGTH=36
###################
# Optional arguments
###################
# WIN=75
# BIN=20
# STAMPIPES=~/stampipes # Location of stampipes installation
###################
# REQUIRED MODULES
###################
# module load samtools
# module load python
# module load bedops
# module load bedtools
###################

FAI ?= $(BWAINDEX).fai

# Window
WIN=75
# Bin interval
BINI=20

BAMFILE ?= $(SAMPLE_NAME).$(GENOME).uniques.sorted.bam
RAW_DENSITY_FILE ?= $(SAMPLE_NAME).$(WIN)_$(BINI).$(GENOME).uniques-density.bed.starch
STAMPIPES ?= ~/stampipes

TMPDIR ?= $(shell pwd)

STARCH_OUT=$(SAMPLE_NAME).$(WIN)_$(BINI).normalized.$(GENOME).uniques-density.bed.starch
BW_OUT=$(SAMPLE_NAME).$(WIN)_$(BINI).normalized.$(GENOME).bw

SCALE=1000000

BED_TMP=$(TMPDIR)/$(SAMPLE_NAME).norm.tmp.bed
WIG_TMP=$(TMPDIR)/$(SAMPLE_NAME).norm.tmp.wig

#ALL_TAGCOUNTS := $(shell samtools view -c ${LIBRARY_NAME}.${GENOME}.uniques.sorted.bam)`
#EXTRANUCLEAR_TAGCOUNTS := $(shell samtools view -c ${LIBRARY_NAME}.${GENOME}.uniques.sorted.bam chrM chrC)

all : $(STARCH_OUT) $(BW_OUT)

$(STARCH_OUT) $(BED_TMP) : $(RAW_DENSITY_FILE) $(BAMFILE)
	unstarch $(RAW_DENSITY_FILE) | \
	  awk -v allcounts=`samtools view -c $(BAMFILE)` -v extranuclear_counts=`samtools view -c $(BAMFILE) chrM chrC` -v scale=$(SCALE) 'BEGIN{ tagcount=allcounts-extranuclear_counts }{z=$$5; n=(z/tagcount)*scale; print $$1 "\t" $$2 "\t" $$3 "\t" $$4 "\t" n }' \
	  | tee $(BED_TMP) \
	  | starch - > $(STARCH_OUT)

$(WIG_TMP) : $(BED_TMP)
	cat $(BED_TMP) \
	  | awk -v s=0 '{ if( $$5 != s ){ print $$0 } }' \
	  | awk 'BEGIN{ OFS = "\t"; chr = ""}{ if( chr == "" ){ chr=$$1; print "variableStep chrom="chr" span=20" } if( $$1 == chr ){ print $$2, $$5 } else { chr=$$1; print "variableStep chrom="chr" span=20"; print $$2, $$5} }' \
	  > $(WIG_TMP)

$(BW_OUT) : $(WIG_TMP)
	wigToBigWig -clip $(WIG_TMP) $(FAI) $(BW_OUT)

