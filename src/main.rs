use clap::Parser;

use hermes_rs::run::cli::Cli;
use hermes_rs::run::runner;

fn main() {
    let cli = Cli::parse();

    if let Err(e) = runner::run(&cli) {
        eprintln!("Error: {e}");
        std::process::exit(1);
    }
}
