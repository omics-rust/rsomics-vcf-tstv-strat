mod cli;
use clap::Parser;
use cli::Cli;
use rsomics_common::Tool;
use std::process::ExitCode;

fn main() -> ExitCode {
    let args = Cli::parse();
    args.run()
}
