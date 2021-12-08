//use rust_htslib::{bam, bam::Read, tpool};
//use fastq::{parse_path, thread_reader, Parser};
//use simple_error::{simple_error, SimpleError};
//use std::collections::HashMap;
//use std::rc::Rc;



//use std::cell::{RefCell};

//use std::error::Error;
//use std::io::Read;
use std::path::{PathBuf};
//use std::str::from_utf8;
use structopt::StructOpt;

//use autocompress::{create, iothread::IoThread, open, CompressionLevel};

//use needletail::parser::parse_fastx_reader;
//use needletail::parser::{parse_fastx_file, parse_fastx_stdin, SequenceRecord};

use demux::{Demultiplexer,BarcodeAssignment};
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
}


fn run(options: &Opt) {
    let demuxer = Demultiplexer {
        assignments: options.barcodes.clone(),
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
