use needletail::parser::{SequenceRecord};
use std::fs::File;
use flate2::write::{GzEncoder};
use flate2::GzBuilder;
use flate2::read::GzDecoder;

use std::path::{PathBuf};
use flate2::Compression;
use std::io::{Read,BufReader};
use autocompress::{Decoder,open_or_stdin};
use autocompress::iothread::{IoThread,ThreadReader,ThreadWriter};

use seq_io;


pub struct IoPool {
    pool: IoThread,
    compression_level: autocompress::CompressionLevel,
}
impl IoPool {
    pub fn new(threads: usize) -> IoPool {
        IoPool { pool: IoThread::new(threads), compression_level: autocompress::CompressionLevel::Default }
    }

    pub fn get_writers(&mut self, filenames: Vec<PathBuf>) -> Vec<ThreadWriter<impl std::io::Write>> {
        let p = &self.pool;
        filenames.iter().map (|filename|{
            p.add_writer(autocompress::create(filename, self.compression_level).expect("Could not open file for writing"))
        }).collect()

    }
    pub fn get_writer(&mut self, filename: PathBuf) -> ThreadWriter<impl std::io::Write> {
        let w = autocompress::create(filename, self.compression_level).expect("Could not open file for writing");
        self.pool.add_writer(w)
    }
    pub fn get_reader(&mut self, filename: Option<PathBuf>) -> seq_io::fastq::Reader<Box<dyn Read + '_>> {

        seq_io::fastq::Reader::new(
            match filename {
                None => Box::new(open_or_stdin(filename).expect("Couldn't open input")),
                Some(x) => Box::new(self.get_raw_reader(x)),

            })
    }
    fn get_raw_reader(&mut self, filename: PathBuf) -> ThreadReader< impl std::io::Read > {
        self.pool.open(filename).expect("Could not open input file")
    }
    // TODO: Figure out how to implement this (err: Read is not Send)
    //fn get_raw_stdin_reader(&mut self) -> ThreadReader< impl std::io::Read > {
    //    self.pool.add_reader(open_or_stdin(None).expect("Couldn't open input")).expect("Couldn't start read thread")
    //}
}

//fn get_raw_file_reader(filename: PathBuf) -> ThreadReader<'_, impl std::io::Read> {
//    let thread_pool = IoThread::new(1);
//    thread_pool.open(filename).expect("Couldn't open file")
//}
//fn get_stdin_reader() {
//    let thread_pool = IoThread::new(1);
//    thread_pool.open(std::io::stdin())
//}

pub struct FastqRawWriter {
    writer: File
}
impl FastqRawWriter {
    pub fn new(path: PathBuf) -> FastqRawWriter {
        FastqRawWriter {
            writer: File::create(path).unwrap(),
        }
    }
    pub fn write(&mut self, record: seq_io::fastq::RefRecord) -> () {
        record.write_unchanged(&mut self.writer).expect("Couldn't write :(")
    }

    pub fn close(&mut self) -> Result<(), std::io::Error> {
        self.writer.sync_all()
    }
}
//impl Write for FastqRawWriter {
//}

pub trait Write {
    fn write(&mut self, record: SequenceRecord) -> ();
    fn close(&mut self) -> Result<(),std::io::Error>;
}

struct FastqGzWriter {
    writer: GzEncoder<File>
}
impl FastqGzWriter {
    pub fn new (path: PathBuf) -> FastqGzWriter {
        FastqGzWriter {
            writer: GzBuilder::new().write(File::create(path).expect("Could not open file for writing"), Compression::default()),
        }
    }
}
impl Write for FastqGzWriter {
    fn write(&mut self, record: SequenceRecord) -> () {
        record.write(&mut self.writer, None).expect("Couldn't write :(")
    }

    fn close(&mut self) -> Result<(), std::io::Error> {
        self.writer.try_finish()
    }
}

pub struct FastqWriter {
    file: File,
}
impl FastqWriter {
    pub fn new(path: PathBuf) -> FastqWriter {
        FastqWriter {
            file: File::create(path).unwrap(),
        }
    }
}
impl Write for FastqWriter {
    fn write(&mut self, record: SequenceRecord) -> () {
        record.write(&mut self.file, None).expect("Couldn't write :(")
    }

    fn close(&mut self) -> Result<(), std::io::Error> {
        self.file.sync_all()
    }
}
