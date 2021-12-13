use flate2::read::GzDecoder;
use std::fs::File;

use autocompress::iothread::{BlockCompress, BlockWriter, IoThread, ThreadReader, ThreadWriter};
use autocompress::{open_or_stdin, CompressionLevel};
use std::io::Read;
use std::path::PathBuf;

use thread_io;

const BUF_SIZE: usize = 256 * 1024;
const QUEUE_LEN: usize = 5;

pub fn with_thread_reader_from_file<F, O, E>(filename: &PathBuf, func: F) -> Result<O, E>
where
    F: FnOnce(&mut thread_io::read::Reader) -> Result<O, E>,
    E: Send,
{
    let f = File::open(filename).unwrap();
    let gz = GzDecoder::new(f);
    //let r = File::open(filename).expect("Could not open input file");
    thread_io::read::reader(BUF_SIZE, QUEUE_LEN, gz, func)
}

pub struct IoPool {
    pool: IoThread,
    compression_level: autocompress::CompressionLevel,
}
impl IoPool {
    pub fn new(threads: usize) -> IoPool {
        IoPool {
            pool: IoThread::new(threads),
            compression_level: autocompress::CompressionLevel::Default,
        }
    }

    pub fn get_writers(
        &mut self,
        filenames: Vec<PathBuf>,
    ) -> Vec<BlockWriter<impl std::io::Write>> {
        let p = &self.pool;
        filenames
            .iter()
            .map(|filename| {
                let gzip_block_compress = BlockCompress::gzip(CompressionLevel::Default);
                let file_writer =
                    std::fs::File::create(filename).expect("Could not open file for writing");

                p.add_block_writer(file_writer, gzip_block_compress)
            })
            .collect()
    }

    pub fn get_writer(&mut self, filename: PathBuf) -> ThreadWriter<impl std::io::Write> {
        let w = autocompress::create(filename, self.compression_level)
            .expect("Could not open file for writing");
        self.pool.add_writer(w)
    }
    pub fn get_reader(
        &mut self,
        filename: Option<PathBuf>,
    ) -> seq_io::fastq::Reader<Box<dyn Read + '_>> {
        seq_io::fastq::Reader::new(match filename {
            None => Box::new(open_or_stdin(filename).expect("Couldn't open input")),
            Some(x) => Box::new(self.get_raw_reader(x)),
        })
    }
    fn get_raw_reader(&mut self, filename: PathBuf) -> ThreadReader<impl std::io::Read> {
        self.pool.open(filename).expect("Could not open input file")
    }
    // TODO: Figure out how to implement this (err: Read is not Send)
    //fn get_raw_stdin_reader(&mut self) -> ThreadReader< impl std::io::Read > {
    //    self.pool.add_reader(open_or_stdin(None).expect("Couldn't open input")).expect("Couldn't start read thread")
    //}
}
