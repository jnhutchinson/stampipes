use crate::config;

use rust_htslib::{bam, bam::Read, tpool};

/// Run the program, all configuration specified in options
pub fn run(options: config::Config) -> Result<(), Box<dyn std::error::Error>> {
    // Setup
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

    // Do the work
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
pub fn move_tag(record: &mut bam::Record, tagname: &[u8], separator: u8) {
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }

    #[test]
    fn move_tag_smoke_test() {
        let tagname = b"RX";
        let sep = '#' as u8;
        let mut basicrecord = bam::Record::new();
        basicrecord.set(b"foo#UMT", None, b"AAAA", &[255, 255, 255, 255]);

        move_tag(&mut basicrecord, tagname, sep);

        assert_eq!(basicrecord.qname(), b"foo");
        match basicrecord.aux(tagname) {
            Some(value) => {
                if let bam::record::Aux::String(v) = value {
                    assert_eq!(v, b"UMT");
                }
            }
            None => {
                panic!("Aux field doesn't exist");
            }
        }
    }
}
