mod config;
mod core;

use crate::config::Config;
use crate::core::run;
use structopt::StructOpt;

fn main() {
    let mut opt = Config::from_args();
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
