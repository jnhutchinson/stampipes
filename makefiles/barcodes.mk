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
	
tmp.barcode_report.txt : $(barcodes)
	$(REPORT_SCRIPT) > $@

barcodes : $(barcodes)
	

%.barcodes.txt : %.fastq.gz
	python $(TALLY_SCRIPT) $^ > $@
