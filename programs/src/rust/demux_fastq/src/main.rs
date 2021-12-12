use std::io::BufRead;

use std::path::PathBuf;
use structopt::StructOpt;

use demux::{BarcodeAssignment, Demultiplexer};
mod demux;
mod fastq;

#[derive(Debug, StructOpt)]
#[structopt(
    name = "demux_fastq",
    about = "Demultiplexes a FastQ file based on barcode"
)]
struct Opt {
    /// Set the number of threads used
    #[structopt(short = "@", long = "threads", default_value = "2")]
    threads: usize,

    /// Set the compression level of the resulting Fastq files.
    #[structopt(short = "z", long = "compression", default_value = "6")]
    compression_level: u32,

    /// Path to the input fastq file.
    /// If omitted or `-`, will read from standard input.
    #[structopt(parse(from_os_str))]
    input: Option<PathBuf>,

    /// Barcodes and corresponding files to write to
    #[structopt(short = "b", long = "barcode")]
    barcodes: Vec<BarcodeAssignment>,

    #[structopt(short = "B", long = "barcode-config")]
    barcode_config_file: Option<PathBuf>,
}

fn parse_barcode_config(reader: impl std::io::Read) -> Vec<BarcodeAssignment> {
    std::io::BufReader::new(reader)
        .lines()
        .map(|l| l.unwrap().parse().unwrap())
        .collect()
}

fn run(options: &Opt) {
    let mut assignments = Vec::new();
    if let Some(f) = options.barcode_config_file.clone() {
        assignments = parse_barcode_config(std::fs::File::open(f).expect("Couldn't open"))
    }
    assignments.append(&mut options.barcodes.clone());
    let demuxer = Demultiplexer {
        assignments: assignments,
        threads: options.threads,
    };

    demuxer.run(&(options.input))
}

fn main() {
    let mut opt = Opt::from_args();
    if opt.compression_level > 9 {
        opt.compression_level = 9;
    }
    if let Some(ref x) = opt.input {
        if x.display().to_string() == "-" {
            opt.input = None
        }
    }
    run(&opt)
}
