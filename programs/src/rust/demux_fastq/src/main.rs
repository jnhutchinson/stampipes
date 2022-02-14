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

    /// Set the number of allowable mismatches (for the whole barcode)
    #[structopt(short = "m", long = "mismatches", default_value = "0")]
    mismatches: u8,

    /// Path to the input fastq file
    /// If omitted or `-`, will read from standard input.
    #[structopt(parse(from_os_str))]
    input: Option<PathBuf>,

    /// Barcodes and corresponding files to write to
    #[structopt(short = "b", long = "barcode")]
    barcodes: Vec<BarcodeAssignment>,

    /// A file containing barcodes and files, one per line. like 'ACTG=out.fq.gz'
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
        mismatches: options.mismatches,
        threads: options.threads,
        assignments,
    };

    demuxer.run(&(options.input))
}

use std::collections::HashMap;
use std::convert::TryInto;

#[derive(Debug)]
struct FieldDecodeError;
fn find_blocksize_in_gzip_header_extra<'a>(extra: &'a [u8]) -> Result<u16, FieldDecodeError> {
    //let total_len = extra.len();
    let mut pos: usize = 0;
    loop {
        match extra[pos..] {
            [b'B', b'C', 2, 0, a1, a2, ..] => return Ok(1 + u16::from_le_bytes([a1, a2])),
            [] => return Err(FieldDecodeError),
            [_, _, s1, s2, ..] => {
                let info_size = u16::from_le_bytes([s1, s2]) as usize;
                pos += info_size + 4;
            }
            _ => return Err(FieldDecodeError),
        };
    }
}
use flate2::GzHeader;
fn get_block_offset(header: GzHeader) {}

use flate2::read::{GzDecoder, MultiGzDecoder};
use std::collections::VecDeque;
use std::fs::{metadata, File};
use std::io::prelude::*;
use std::io::Read;
use std::io::SeekFrom;
struct BGZFReader {
    f: File,
    threads: usize,
    filename: PathBuf,
    block_locations: VecDeque<(u64, u64)>,
    readers: VecDeque<DecompressedBlock>,
    //cur_reader: Option<BufReader<GzDecoder<File>>>
}

use std::io::BufReader;
use std::sync::Arc;
use std::sync::Mutex;

struct DecompressedBlock {
    buf: Arc<Mutex<BufReader<MultiGzDecoder<File>>>>,
}

fn merge_blocks(blocks: VecDeque<(u64, u64)>, chunk_size: u64) -> VecDeque<(u64, u64)> {
    if blocks.len() <= 1 || chunk_size <= 1 {
        return blocks;
    }
    println!("in: {:?}", blocks);
    let mut out = VecDeque::new();
    let mut cur_block = blocks[0];
    let mut count = 1;

    for block in blocks.into_iter().skip(1) {
        count += 1;
        if count <= chunk_size {
            cur_block.1 += block.1;
        } else {
            out.push_back(cur_block);
            count = 1;
            cur_block = block;
        }
    }

    //if count > 1 {
    out.push_back(cur_block);
    //}
    println!("out: {:?}", out);
    out
}

const MERGE_LEVEL: usize = 2;

use std::thread;
impl BGZFReader {
    pub fn new(filename: &std::path::Path) -> BGZFReader {
        let f = File::open(filename.clone()).unwrap();
        let mut offset: u64 = 0;
        let metadata = f.metadata().unwrap();
        let file_len = metadata.len();

        let mut block_locations = VecDeque::new();
        while offset < file_len {
            let mut f = File::open(filename.clone()).unwrap();
            f.seek(SeekFrom::Start(offset)).unwrap();
            let gz = GzDecoder::new(f);
            let h = gz.header().unwrap();
            let blocksize = find_blocksize_in_gzip_header_extra(h.extra().unwrap()).unwrap();
            block_locations.push_back((offset, blocksize as u64));
            offset += blocksize as u64;
        }
        let readers = VecDeque::with_capacity(8);
        //let new_locs =
        BGZFReader {
            f: File::open(filename).unwrap(),
            filename: filename.into(),
            //block_locations: block_locations,
            block_locations: merge_blocks(block_locations, MERGE_LEVEL as u64),
            readers,
            threads: 8,
            //cur_reader: None,
        }
    }

    fn start_reader(&mut self) {
        match self.block_locations.pop_front() {
            None => (),
            Some((offset, _size)) => {
                let mut f = File::open(self.filename.clone()).unwrap();
                f.seek(SeekFrom::Start(offset)).unwrap();
                let gz = MultiGzDecoder::new(f);
                let br = BufReader::with_capacity(MERGE_LEVEL * 65536, gz);
                let buffer = Arc::new(Mutex::new(br));
                //m.lock();
                let thread_buf = Arc::clone(&buffer);
                self.readers.push_back(DecompressedBlock { buf: buffer });
                // How do we not spawn a new thread every time?
                // TODO: Worker pool!
                thread::spawn(move || {
                    thread_buf.lock().unwrap().fill_buf().unwrap();
                });
            }
        }
    }
    //pub fn get_reader(self) -> Box<impl Read> {
    //Box::new(self)
    //}
}
use std::sync::MutexGuard;
impl Read for BGZFReader {
    fn read(&mut self, into: &mut [u8]) -> std::io::Result<usize> {
        // Spin up threads to read
        while self.readers.len() < self.threads && !self.block_locations.is_empty() {
            self.start_reader();
        }

        if self.readers.is_empty() {
            return Ok(0);
        }
        //let cur_reader = self.readers[0].buf;

        //if self.cur_reader.is_none() {
        //    match self.readers.pop_front() {
        //        None => return Ok(0),
        //        Some(r) => {
        //            let buffer = r.buf.lock().unwrap();
        //            self.cur_reader = Some(buffer);
        //        },
        //    }
        //}
        let b = &self.readers[0].buf;
        match b.lock().unwrap().read(into) {
            Ok(0) => {} // Out fallthrough to below
            Ok(n) => return Ok(n),
            Err(e) => return Err(e),
        };
        // This reader empty, yeet!
        self.readers.pop_front();
        self.read(into)
    }
}

fn run_debug(options: &Opt) {
    //let br = BGZFReader::new(&options.input.clone().unwrap());
    use rust_htslib::{bgzf, tpool};
    let pool = tpool::ThreadPool::new(options.threads as u32).unwrap();
    let mut br = bgzf::Reader::from_path(options.input.clone().unwrap()).unwrap();
    br.set_thread_pool(&pool).unwrap();
    //let br = fastq::bgzip_reader(options.input.clone(), Some(pool));

    let mut count = 0;
    for read in seq_io::fastq::Reader::new(br).records() {
        count += 1;
        if count % 100000 == 0 {
            println!("read:  {}", count);
        }
    }
    println!("Total: {}", count);

    //let mut offset: u64 = 0;
    //let file_len = metadata(options.input.clone().unwrap()).unwrap().len();

    //while offset < file_len {
    //    let mut f = File::open(options.input.clone().unwrap()).unwrap();
    //    f.seek(SeekFrom::Start(offset)).unwrap();
    //    let gz = GzDecoder::new(f);
    //    let h = gz.header().unwrap();
    //    let bl2 = find_blocksize_in_gzip_header_extra(h.extra().unwrap()).unwrap();
    //    println!("{:?}", bl2);
    //    offset += bl2 as u64;
    //}

    println!("done!");
    //let gz = GzDecoder::new(f);
}

fn main() {
    let mut opt = Opt::from_args();
    if let Some(ref x) = opt.input {
        if x.display().to_string() == "-" {
            opt.input = None
        }
    }
    // use  rust_htslib::bgzf::Reader;
    // use rust_htslib::tpool::ThreadPool;
    // let mut r = Reader::from_path(opt.input.unwrap().clone()).unwrap();

    // //r.set_thread_pool(&ThreadPool::new(1).unwrap());

    // let mut count = 0;
    // for read in seq_io::fastq::Reader::new(r).records(){
    //     count += 1;
    // }
    // println!("Count: {}", count);

    //run_debug(&opt)
    run(&opt)
    //
}
