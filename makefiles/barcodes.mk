TALLY_SCRIPT ?= $(STAMPIPES)/scripts/flowcells/tallybarcodes.py
REPORT_SCRIPT ?= $(STAMPIPES)/scripts/flowcells/barcode_report.sh

SAMPLE_FASTQ = $(shell ls Project_*/Sample*/*_R1_???.fastq.gz)
UNDET_FASTQ = $(shell ls Undetermined_indices/Sample*/*_R1_???.fastq.gz)
FASTQ ?= $(SAMPLE_FASTQ) $(UNDET_FASTQ)
barcodes = $(patsubst %.fastq.gz,%.barcodes.txt,$(FASTQ))

all : info barcodes report

info :
	@echo "------"
	@echo "PROGRAM VERSIONS"
	@echo "------"
	@echo "METADATA"
	@echo "------"
	@echo "SAMPLE_NAME: " $(SAMPLE_NAME)
	@echo "------"
	@echo "FASTQ_COUNT: " $(words $(FASTQ))
	@echo "------"

report: barcode_report.txt
	
barcode_report.txt : $(barcodes)
	SGE_RREQ=" -N .tb$(FLOWCELL)-report " $(REPORT_SCRIPT) > $@

barcodes : $(barcodes)
	

%.barcodes.txt : %.fastq.gz
	SGE_RREQ=" -N .tb$(FLOWCELL)-fasta " python $(TALLY_SCRIPT) $^ > $@
