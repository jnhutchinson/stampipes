use std::path::PathBuf;
use structopt::StructOpt;

/// Opt contains all of our configuration for this program
#[derive(Debug, StructOpt)]
#[structopt(name = "umt_tagger", about = "Moves UMTs from the read name to a tag")]
pub struct Config {
    /// Set the number of threads used for compression.
    ///
    /// These threads are shared between the reader and writer.
    /// There is one additional "main" thread used.
    #[structopt(short = "@", long = "threads", default_value = "2")]
    pub threads: u32,

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
    pub compression_level: u32,

    /// Path to the input BAM file.
    /// If omitted or `-`, will read from standard input.
    #[structopt(parse(from_os_str))]
    pub input: Option<PathBuf>,

    /// Path to the output BAM file to be created.
    /// If omitted or `-`, will read from standard input.
    #[structopt(parse(from_os_str))]
    pub output: Option<PathBuf>,

    /// Read separator character
    #[structopt(short = "s", long = "separator", default_value = "#")]
    pub separator: char,

    /// Tagname to use
    #[structopt(short = "t", long = "tagname", default_value = "RX")]
    pub tagname: String,
}
