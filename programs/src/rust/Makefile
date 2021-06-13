all: release linux

release:
	cargo build --release

debug:
	cargo build

linux:
	TARGET_CC=x86_64-linux-musl-gcc \
	cargo build --release --target x86_64-unknown-linux-musl


