EXE_NAMES = demux_fastq
BINDIR = bin
BINARIES = $(addprefix $(BINDIR)/,$(EXE_NAMES))


.PHONY: all

all: $(BINARIES)


# Go build files
$(BINDIR)/%: src/go/%/*.go
	cd "src/go/$*" && \
	go build . && \
	cd ../../.. && \
	mv "src/go/$*/$*" $@
