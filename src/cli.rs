use clap::{Parser, ValueEnum};
use rsomics_common::{CommonFlags, Result, RsomicsError, Tool, ToolMeta};
use rsomics_vcf_tstv_strat::{tstv_by_count, tstv_by_qual, write_count_table, write_qual_table};
use std::io::{self, BufReader, Write};
use std::path::PathBuf;

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum StratMode {
    /// Stratify Ts/Tv by ALT allele count across samples.
    Count,
    /// Stratify Ts/Tv by QUAL threshold.
    Qual,
}

#[derive(Parser, Debug)]
#[command(
    name = "rsomics-vcf-tstv-strat",
    version,
    about = "Ts/Tv stratified by ALT allele count or QUAL (vcftools --TsTv-by-count / --TsTv-by-qual)",
    long_about = None,
    disable_help_flag = true
)]
pub struct Cli {
    /// Stratification mode: count (by ALT allele count) or qual (by QUAL threshold).
    #[arg(long, short = 'b', value_name = "MODE", default_value = "count")]
    pub by: StratMode,

    /// Input VCF file (omit or use - for stdin).
    #[arg(value_name = "INPUT")]
    pub input: Option<PathBuf>,

    #[command(flatten)]
    pub common: CommonFlags,
}

impl Cli {
    pub fn execute(self) -> Result<()> {
        let stdout = io::stdout();
        let mut out = stdout.lock();

        macro_rules! open_reader {
            ($path:expr) => {
                BufReader::new(
                    std::fs::File::open($path).map_err(|e| {
                        RsomicsError::InvalidInput(format!("{}: {e}", $path.display()))
                    })?,
                )
            };
        }

        let is_stdin = self.input.as_deref().is_none_or(|p| p.as_os_str() == "-");

        match self.by {
            StratMode::Count => {
                let rows = if is_stdin {
                    tstv_by_count(BufReader::new(io::stdin().lock()))?
                } else {
                    tstv_by_count(open_reader!(self.input.as_ref().unwrap()))?
                };
                if self.common.json {
                    let json = serde_json::to_string_pretty(&rows)
                        .map_err(|e| RsomicsError::InvalidInput(e.to_string()))?;
                    writeln!(out, "{json}").map_err(RsomicsError::Io)?;
                } else {
                    write_count_table(&rows, &mut out).map_err(RsomicsError::Io)?;
                }
            }
            StratMode::Qual => {
                let rows = if is_stdin {
                    tstv_by_qual(BufReader::new(io::stdin().lock()))?
                } else {
                    tstv_by_qual(open_reader!(self.input.as_ref().unwrap()))?
                };
                if self.common.json {
                    let json = serde_json::to_string_pretty(&rows)
                        .map_err(|e| RsomicsError::InvalidInput(e.to_string()))?;
                    writeln!(out, "{json}").map_err(RsomicsError::Io)?;
                } else {
                    write_qual_table(&rows, &mut out).map_err(RsomicsError::Io)?;
                }
            }
        }

        Ok(())
    }
}

impl Tool for Cli {
    fn meta() -> ToolMeta {
        META
    }

    fn common(&self) -> &CommonFlags {
        &self.common
    }

    fn execute(self) -> Result<()> {
        self.execute()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;

    #[test]
    fn cli_debug_assert() {
        Cli::command().debug_assert();
    }
}
