use rust_htslib::{bam, bam::Read, tpool};
use std::path::PathBuf;
use structopt::StructOpt;

/// Opt contains all of our configuration for this program
#[derive(Debug, StructOpt)]
#[structopt(name = "umt_tagger", about = "Moves UMTs from the read name to a tag")]
struct Opt {
    /// Set the number of threads used for compression.
    ///
    /// These threads are shared between the reader and writer.
    /// There is one additional "main" thread used.
    #[structopt(short = "@", long = "threads", default_value = "2")]
    threads: u32,

    #[structopt(verbatim_doc_comment)]
    /// Set the compression level of the resulting BAM file.
    ///
    /// 0: No compression (similar to `samtools view -u`)
    /// 1: Fastest compression
    /// 9: Best compression
    #[structopt(
        short = "z",
        long = "compression",
        value_name = "level",
        default_value = "6"
    )]
    compression_level: u32,

    /// Path to the input BAM file.
    /// If omitted or `-`, will read from standard input.
    #[structopt(parse(from_os_str))]
    input: Option<PathBuf>,

    /// Path to the output BAM file to be created.
    /// If omitted or `-`, will read from standard input.
    #[structopt(parse(from_os_str))]
    output: Option<PathBuf>,

    /// Read separator character
    #[structopt(short = "s", long = "separator", default_value = "#")]
    separator: char,

    /// Tagname to use
    #[structopt(short = "t", long = "tagname", default_value = "RX")]
    tagname: String,
}

/// Run the program, all configuration specified in options
fn run(options: Opt) -> Result<(), Box<dyn std::error::Error>> {
    let mut bam = match options.input {
        Some(path) => bam::Reader::from_path(path)?,
        None => bam::Reader::from_stdin()?,
    };

    let header = bam::Header::from_template(bam.header());
    let mut output = match options.output {
        Some(path) => bam::Writer::from_path(path, &header, bam::Format::BAM)?,
        None => bam::Writer::from_stdout(&header, bam::Format::BAM)?,
    };

    let threadpool = tpool::ThreadPool::new(options.threads)?;
    bam.set_thread_pool(&threadpool)?;
    output.set_thread_pool(&threadpool)?;
    output.set_compression_level(bam::CompressionLevel::Level(options.compression_level))?;

    let tagname = options.tagname.as_bytes();
    let separator = options.separator as u8;

    let mut record = bam::Record::new();
    while let Some(result) = bam.read(&mut record) {
        result.expect("BAM parsing failed");

        move_tag(&mut record, tagname, separator);
        output.write(&record)?;
    }
    Ok(())
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
        let umt = &owned[sep_location + 1..];
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

    let result = run(opt);

    std::process::exit(match result {
        Ok(_) => 0,
        Err(err) => {
            eprintln!("error: {:}", err);
            1
        }
    });
}
