//use rust_htslib::{bam, bam::Read, tpool};
//use fastq::{parse_path, thread_reader, Parser};
use simple_error::{simple_error, SimpleError};
use std::collections::HashMap;
use std::fs::File;
use std::rc::Rc;

use std::cell::{RefCell};

//use std::error::Error;
//use std::io::Read;
use std::path::{PathBuf};
//use std::str::from_utf8;
use structopt::StructOpt;

//use autocompress::{create, iothread::IoThread, open, CompressionLevel};

//use needletail::parser::parse_fastx_reader;
use needletail::parser::{parse_fastx_file, parse_fastx_stdin, SequenceRecord};

pub struct Demultiplexer {
    assignments: Vec<BarcodeAssignment>,
}
impl Demultiplexer {
    fn run(&self, input: &Option<PathBuf>) {
        let (demux_map, mut writers) = generate_demux_map(&self.assignments, 0);

        let mut reader = match input {
            Some(filename) => {
                parse_fastx_file(filename).expect("couldn't open input file as fastq")
            }
            None => parse_fastx_stdin().expect("invalid reader"),
        };

        while let Some(r) = reader.next() {
            let record = r.expect("Invalid record");
            let id = record.id();
            let sep = ':' as u8;
            let slicestart = id.len()
                - id.iter()
                    .rev()
                    .position(|x| *x == sep)
                    .expect("no bc in record");
            let bc = &id[slicestart..];

            match demux_map.get(bc) {
                Some(index) => record.write(&mut writers[*index], None).unwrap(),
                None => (),
            };
        }
    }
}

#[derive(Debug, Clone)]
struct BarcodeAssignment {
    barcode: String,
    filepath: PathBuf,
}
impl std::str::FromStr for BarcodeAssignment {
    type Err = SimpleError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.split_once('=') {
            Some((barcode, path_string)) if !barcode.is_empty() && !path_string.is_empty() => {
                Ok(BarcodeAssignment {
                    barcode: barcode.to_string(),
                    filepath: PathBuf::from_str(path_string).unwrap(),
                })
            }

            _ => Err(simple_error!(
                "'{}' does not match barcode format: 'barcode=file'",
                s
            )),
        }
    }
}

//#[derive(Clone)]
struct FastqWriter {
    path: PathBuf,
    file: File,
}
impl FastqWriter {
    fn new(path: PathBuf) -> FastqWriter {
        FastqWriter {
            path: path.clone(),
            file: File::create(path).unwrap(),
        }
    }

    fn write(&mut self, record: SequenceRecord) -> () {
        record.write(&mut self.file, None).unwrap()
    }
}

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

/// This returns a map of barcodes to Open writers
fn generate_demux_map(
    barcodes: &Vec<BarcodeAssignment>,
    mismatches: u8,
) -> (HashMap<Vec<u8>, usize>, Vec<File>) {
    let mut path_map: HashMap<PathBuf, usize> = HashMap::new();

    let mut v = Vec::new();
    {
        let mut i = 0;
        for ba in barcodes.iter() {
            let path = ba.filepath.clone();
            if !path_map.contains_key(&path) {
                let writer =
                    File::create(ba.filepath.clone()).expect("Couldn't open file for writing");
                v.push(writer);
                path_map.insert(path, i);
                i += 1;
            }
        }
    }

    let mut barcode_map = HashMap::new();
    for ba in barcodes.iter() {
        let w = path_map[&ba.filepath];
        let b = ba.barcode.as_bytes().to_vec();
        barcode_map.insert(b, w);
    }
    (barcode_map, v)
}

fn run(options: &Opt) {
    let demuxer = Demultiplexer {
        assignments: options.barcodes.clone(),
    };

    demuxer.run(&(options.input))
    /*
    parse_path(options.input.as_ref(), |parser| {
        let results: Vec<usize> = parser
            .parallel_each(options.threads, |record_sets| {
                let mut count = 0;
                for record_set in record_sets {
                    for record in record_set.iter() {
                        count += 1;
                        let head = fastq::Record::head(&record);
                        println!("{}", String::from_utf8_lossy(head));
                    }
                }
                count
            })
            .expect("Ah heck some error");
        //println!("{}", results.iter().sum::<usize>());
    })
    .expect("Invalid compression");
    */
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
