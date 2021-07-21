all: release linux

.PHONY: all check clean debug linux releae test

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
