all: release linux
.PHONY: all check clean debug linux release test pgo

PGO_DIR=$(shell pwd)/tmp/pgo

debug:
	cargo build

release:
	cargo build --release

linux:
	TARGET_CC=x86_64-linux-musl-gcc \
	cargo build --release --target x86_64-unknown-linux-musl

test:
	cargo test

clean:
	rm -rf target/

check:
	cargo check

pgo:
	# Step 0 - remove old profiling data
	rm -rf "$(PGO_DIR)"
	mkdir -p "$(PGO_DIR)"
	# Step 1 - compile an instrumented version of the binary
	RUSTFLAGS="-Cprofile-generate=$(PGO_DIR)" \
	    cargo build --release # --target=x86_64-unknown-linux-gnu
	# Step 2: Run several times
	./target/release/demux_fastq test_data/lorg.fq.bgz --barcode-config test_data/barcode_config.2 --threads 2
	./target/release/demux_fastq test_data/lorg.fq.bgz --barcode-config test_data/barcode_config.10 --threads 5
	./target/release/demux_fastq test_data/lorg.fq.bgz --barcode-config test_data/barcode_config.all --threads 10
	#./target/x86_64-unknown-linux-gnu/release/demux_fastq
	#./target/x86_64-unknown-linux-gnu/release/demux_fastq
	# Step 3: Merge profile data
	rust-profdata merge -o "$(PGO_DIR)/merged.profdata" "$(PGO_DIR)"
	# Step 4: Use the `.profdata` file for guiding optimizations
	RUSTFLAGS="-Cprofile-use=$(PGO_DIR)/merged.profdata" \
	    cargo build --release
