###############
# These variables should be passed in for the makefile to work properly
###############
# INPUT_BAM_FILE=
# OUTPUT_BAM_FILE=
###############

TMPDIR ?= $(shell pwd)

NSLOTS ?= 1
THREADS ?= $(NSLOTS)

DEFAULT_MEM ?= 4
MAX_MEM ?= 4

SAMTOOLS_FILTER_OPTIONS ?= -f 2 -F 512

# Check if we can sort in parallel

SORT := sort

SORT_PARALLEL = $(shell echo `sort --version | grep ^sort | sed 's/^.* //g'` \>= 8.5 | bc)
ifeq ($(SORT_PARALLEL), 1)
	SORT += --parallel=$(THREADS)
endif

all: markdups

markdups: $(OUTPUT_BAM_FILE)

$(OUTPUT_BAM_FILE): $(TMPDIR)/reads.markdups.sam
	samtools view -Sbu $^ > $@

$(TMPDIR)/reads.markdups.sam: $(TMPDIR)/reads.sam $(TMPDIR)/reads.dups
	samtools view -SH $(TMPDIR)/reads.sam > $@
	
	samtools view -S $(TMPDIR)/reads.sam \
		| join -j 1 -v 1 - $(TMPDIR)/reads.dups \
		| tr " " "\t" \
		| awk -v OFS="\t" '{ $$2 = and($$2, compl(lshift(1, 10))); print; }' \
	>> $@

	samtools view -S $(TMPDIR)/reads.sam \
		| join -j 1 - $(TMPDIR)/reads.dups \
		| tr " " "\t" \
		| awk -v OFS="\t" '{ $$2 = or($$2, lshift(1, 10)); print; }' \
	>> $@

# Prints out all duplicates (except the first one) (previous step
# sorts reads by highest MAPQ
#
# Strategy is to find all the duplicate lines, then remove the first
# instance of a duplicate, and return the rest

$(TMPDIR)/reads.dups: $(TMPDIR)/reads.alldups $(TMPDIR)/reads.firstdup
	join -j 1 -v 1 $(TMPDIR)/reads.alldups $(TMPDIR)/reads.firstdup > $@ 

$(TMPDIR)/reads.firstdup: $(TMPDIR)/reads.sorted
	cat $^ | uniq -f2 -d | cut -f1 | $(SORT) --buffer-size=$(MAX_MEM)G > $@

$(TMPDIR)/reads.alldups: $(TMPDIR)/reads.sorted
	cat $^ | uniq -f2 -D | cut -f1 | $(SORT) --buffer-size=$(MAX_MEM)G > $@

# Sort the fragments by position, strand, UMI and then quality
# TAG ID, MAPQ, chr, start, end, strand, UMI

$(TMPDIR)/reads.sorted: $(TMPDIR)/reads.sam
	samtools view $(SAMTOOLS_FILTER_OPTIONS) -S $^ \
		| awk -f $(STAMPIPES)/scripts/umi/umi_sort_sam_annotate.awk \
		| $(SORT) --buffer-size=$(MAX_MEM)G -k3,3 -k4,4g -k5,5g -k6,6 -k7,7 -k2,2gr \
	>$@

# Sort by name to start

$(TMPDIR)/reads.sam: $(INPUT_BAM_FILE)
	samtools view -H $^ > $@
	samtools view $^ | $(SORT) --buffer-size=$(MAX_MEM)G -k1,1 >> $@

