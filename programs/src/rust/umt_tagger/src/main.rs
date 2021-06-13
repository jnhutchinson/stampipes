use rust_htslib::{bam, bam::Read};
use std::env;

fn main() {
    let args: Vec<String> = env::args().collect();
    let input_name = &args[1];
    let output_name = &args[2];

    let mut bam = bam::Reader::from_path(input_name).unwrap();
    let header = bam::Header::from_template(bam.header());

    let mut output = bam::Writer::from_path(output_name, &header, bam::Format::BAM).unwrap();
    output.set_threads(8).expect("Could not set threads");

    let mut record = bam::Record::new();
    while let Some(result) = bam.read(&mut record) {
        match result {
            Ok(_) => {
                // Get a copy of the qname
                let qname = record.qname();

                let position = qname.iter().position(|x| x == &('#' as u8));
                match position {
                    Some(x) => {
                        let owned = qname.to_owned().into_boxed_slice();
                        let (new_name, mut umt) = owned.split_at(x);
                        // Get rid of delimiter
                        umt = &umt[1..umt.len()];
                        record.set_qname(new_name);
                        record.push_aux("RX".as_bytes(), &bam::record::Aux::String(umt))
                    }
                    None => ()
                }

                output.write(&record).expect("Could not write read");
            },
            Err(_) => panic!("BAM parsing failed...")
        }
    }
}
