use std::path::PathBuf;
use rust_htslib::{bam, bam::Read, tpool};
use structopt::StructOpt;

/// Opt contains all of our configuration for this program
#[derive(Debug, StructOpt)]
#[structopt(name = "umt_tagger", about = "Moves UMTs from the read name to a tag")]
struct Opt {
    /// Set the number of threads used by the underlying htslib.
    /// These threads are shared between the reader and writer.
    /// There is one additional "main" thread used to do the work of moving the UMTs.
    #[structopt(short = "@", long = "threads", default_value = "2")]
    threads: u32,

    /// Set the compression level of the resulting BAM file.
    /// 0: No compression (similar to `samtools view -u`)
    /// 1: Fastest compression
    /// 9: Best compression
    #[structopt(short = "z", long = "compression", default_value = "6")]
    compression_level: u32,

    /// Path to the input BAM file.
    /// If omitted or `-`, will read from standard input.
    #[structopt(parse(from_os_str))]
    input: Option<PathBuf>,

    /// Path to the output BAM file to be created.
    /// If omitted or `-`, will read from standard input.
    // Output file, stdout if not present
    #[structopt(parse(from_os_str))]
    output: Option<PathBuf>,

    /// Read separator character
    #[structopt(short = "s", long = "separator", default_value = "#")]
    separator: char,

    /// Tagname to use
    #[structopt(short = "tag", long = "tagname", default_value = "RX")]
    tagname: String,
}


/// Run the program, all configuration specified in options
fn run(options: Opt) {
    let mut bam = match options.input {
        Some(path) => bam::Reader::from_path(path).expect("Could not open input file"),
        None => bam::Reader::from_stdin().expect("Could not open stdin"),
    };

    let header = bam::Header::from_template(bam.header());
    let mut output = match options.output {
        Some(path) => bam::Writer::from_path(path, &header, bam::Format::BAM).expect("Could not open output file"),
        None => bam::Writer::from_stdout(&header, bam::Format::BAM).expect("Could not open stdout"),
    };

    let threadpool = tpool::ThreadPool::new(options.threads).expect("Could not create threadpool");
    bam.set_thread_pool(&threadpool)
        .expect("Could not modify input threadpool");
    output.set_thread_pool(&threadpool)
        .expect("Could not modify output threadpool");
    output.set_compression_level(bam::CompressionLevel::Level(options.compression_level))
        .expect("Could not set output compression level");

    let tagname = options.tagname.as_bytes();
    let separator = options.separator as u8;

    let mut record = bam::Record::new();
    while let Some(result) = bam.read(&mut record) {
        match result {
            Ok(_) => {
                move_tag(&mut record, tagname, separator);
                output.write(&record).expect("Could not write read");
            },
            Err(_) => panic!("BAM parsing failed...")
        }
    }
}

/// Checks for the presence of the UMT separator in the read name
/// If found, moves everything past the separator into the specified tag, and shortens up the name.
fn move_tag(record: &mut bam::Record, tagname: &[u8], separator: u8) {
    let qname = record.qname();
    let position_from_end = qname.iter().rev().position(|x| *x == separator);
    if let Some(x) = position_from_end {
        // Get a copy of the qname that we can edit
        let owned = qname.to_owned().into_boxed_slice();
        let sep_location = owned.len() - x - 1;

        let new_name = &owned[..sep_location];
        let umt = &owned[sep_location+1..];
        record.set_qname(new_name);
        record.push_aux(tagname, &bam::record::Aux::String(umt))
    }
}

fn main() {
    let mut opt = Opt::from_args();
    // Fix up options...
    if opt.compression_level > 9 {
        opt.compression_level = 9;
    }

    // Need at least one thread
    if opt.threads == 0 {
        opt.threads = 1;
    }

    run(opt);
}
