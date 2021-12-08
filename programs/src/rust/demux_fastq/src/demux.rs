use std::collections::HashMap;
use std::io::BufReader;
use simple_error::{simple_error, SimpleError};
use std::path::{PathBuf};
use std::rc::Rc;


//use crate::fastq::Write;
use seq_io::fastq::{Reader,Record};

use crate::fastq::{FastqRawWriter,IoPool};

use std::io::Write;
use autocompress::{Decoder,open_or_stdin};
use seq_io;
//mod fastq;


pub struct Demultiplexer {
    pub assignments: Vec<BarcodeAssignment>,
}
impl Demultiplexer {
    pub fn run(&self, input: &Option<PathBuf>) {

        let mut iopool = IoPool::new(5);

       // let mut reader = iopool.get_reader(input.clone());
        let mut reader = seq_io::fastq::Reader::new(
            open_or_stdin(input.clone()).expect("Couldn't open input")
        );

        let (demux_map, mut writers) = generate_demux_map(&self.assignments, &mut iopool, 0);


        while let Some(r) = reader.next() {
            let record = r.expect("Invalid record");
            let id = record.head();
            let sep = ':' as u8;
            let slicestart = id.len()
                - id.iter()
                    .rev()
                    .position(|x| *x == sep)
                    .expect("no bc in record");
            let bc = &id[slicestart..];

            match demux_map.get(bc) {
                Some(index) => {
                    let w = &mut writers[*index];
                    record.write_unchanged(w).expect("could not write record to output");
                },
                None => (),
            };
        }
        for mut w in writers {
            w.close().expect("Couldn't write");
        }
    }
}

#[derive(Debug, Clone)]
pub struct BarcodeAssignment {
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
/// This returns a map of barcodes to Open writers
fn generate_demux_map<'a>(
    barcodes: &Vec<BarcodeAssignment>,
    pool: &'a mut IoPool,
    _mismatches: u8,
) -> (HashMap<Vec<u8>, usize>,  Vec<autocompress::iothread::ThreadWriter<'a, impl std::io::Write>>) {
    let mut path_map: HashMap<PathBuf, usize> = HashMap::new();

    let mut paths = Vec::new();
    {
        let mut i = 0;
        for ba in barcodes.iter() {
            let path = ba.filepath.clone();
            if !path_map.contains_key(&path) {
                paths.push(path.clone());
                //let writer =
                    //pool.get_writer(ba.filepath.clone());
                    //FastqRawWriter::new(ba.filepath.clone());
                    //File::create(ba.filepath.clone()).expect("Couldn't open file for writing");
                //v.push(writer);
                path_map.insert(path, i);
                i += 1;
            }
        }
    }
    let v = pool.get_writers(paths);

    let mut barcode_map = HashMap::new();
    for ba in barcodes.iter() {
        let w = path_map[&ba.filepath];
        let b = ba.barcode.as_bytes().to_vec();
        barcode_map.insert(b, w);
    }
    (barcode_map, v)
}
