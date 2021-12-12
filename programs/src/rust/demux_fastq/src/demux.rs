use simple_error::{simple_error, SimpleError};
use std::collections::HashMap;
use std::path::PathBuf;

use seq_io::fastq::{Record, RefRecord};

use crate::fastq::{with_thread_reader_from_file, IoPool};

use seq_io;

pub struct Demultiplexer {
    pub assignments: Vec<BarcodeAssignment>,
    pub threads: usize,
}
impl Demultiplexer {
    pub fn run(&self, input: &Option<PathBuf>) {
        let mut iopool = IoPool::new(self.threads);
        let (demux_map, mut writers) = generate_demux_map(&self.assignments, &mut iopool, 0);
        let result = with_thread_reader_from_file(input.as_ref().unwrap(), |rdr| {
            let mut reader = seq_io::fastq::Reader::new(rdr);
            while let Some(r) = reader.next() {
                let record = r.expect("Invalid record");
                let bc = get_record_bc(&record);
                match demux_map.get(bc) {
                    Some(index) => {
                        let w = &mut writers[*index];
                        record
                            .write_unchanged(w)
                            .expect("could not write record to output");
                    }
                    None => (),
                };
            }
            for mut w in writers {
                w.close().expect("Couldn't write");
            }
            Ok::<_, std::io::Error>(())
        });
        result.unwrap();
    }
}

fn get_record_bc<'a>(record: &'a RefRecord) -> &'a [u8] {
    let id = record.head();
    let sep = ':' as u8;
    let slicestart = id.len()
        - id.iter()
            .rev()
            .position(|x| *x == sep)
            .expect("no bc in record");
    let bc = &id[slicestart..];
    bc
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
) -> (
    HashMap<Vec<u8>, usize>,
    Vec<autocompress::iothread::BlockWriter<'a, impl std::io::Write>>,
) {
    let mut path_map: HashMap<PathBuf, usize> = HashMap::new();

    let mut paths = Vec::new();
    {
        let mut i = 0;
        for ba in barcodes.iter() {
            let path = ba.filepath.clone();
            if !path_map.contains_key(&path) {
                paths.push(path.clone());
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
