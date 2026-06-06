// pattern: Functional Core

use std::error::Error;
use std::fmt::{self, Display};
use std::path::PathBuf;

/// Error type for genotype source IO, representation, filtering, and matrix contract failures.
#[derive(Debug)]
pub enum GenoioError {
    Io {
        path: PathBuf,
        source: std::io::Error,
    },
    InvalidSource {
        path: PathBuf,
        message: String,
    },
    UnsupportedRepresentation {
        message: String,
    },
    SampleFilter {
        requested: usize,
        retained: usize,
        missing: usize,
    },
    MissingData {
        message: String,
    },
    InvalidFilter {
        message: String,
    },
    InternalContract {
        message: String,
    },
}

impl GenoioError {
    /// Build an invalid-source error for a path or logical component.
    pub fn invalid_source(path: impl Into<PathBuf>, message: impl Into<String>) -> Self {
        Self::InvalidSource {
            path: path.into(),
            message: message.into(),
        }
    }

    /// Build an unsupported-representation error.
    pub fn unsupported(message: impl Into<String>) -> Self {
        Self::UnsupportedRepresentation {
            message: message.into(),
        }
    }

    /// Build a sample-filter error with structured counts.
    pub fn sample_filter(requested: usize, retained: usize, missing: usize) -> Self {
        Self::SampleFilter {
            requested,
            retained,
            missing,
        }
    }

    /// Build a missing-data error.
    pub fn missing_data(message: impl Into<String>) -> Self {
        Self::MissingData {
            message: message.into(),
        }
    }

    /// Build an invalid-filter error.
    pub fn invalid_filter(message: impl Into<String>) -> Self {
        Self::InvalidFilter {
            message: message.into(),
        }
    }

    /// Build an internal contract error for impossible states at public boundaries.
    pub fn internal_contract(message: impl Into<String>) -> Self {
        Self::InternalContract {
            message: message.into(),
        }
    }
}

impl Display for GenoioError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io { path, source } => {
                write!(
                    formatter,
                    "failed to read source {}: {source}",
                    path.display()
                )
            }
            Self::InvalidSource { path, message } => {
                write!(
                    formatter,
                    "invalid source {}: {message}",
                    path.display()
                )
            }
            Self::UnsupportedRepresentation { message } => formatter.write_str(message),
            Self::SampleFilter {
                requested,
                retained,
                missing,
            } => write!(
                formatter,
                "missing requested sample(s): requested={requested} retained={retained} missing={missing}"
            ),
            Self::MissingData { message }
            | Self::InvalidFilter { message }
            | Self::InternalContract { message } => formatter.write_str(message),
        }
    }
}

impl Error for GenoioError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            Self::Io { source, .. } => Some(source),
            Self::InvalidSource { .. }
            | Self::UnsupportedRepresentation { .. }
            | Self::SampleFilter { .. }
            | Self::MissingData { .. }
            | Self::InvalidFilter { .. }
            | Self::InternalContract { .. } => None,
        }
    }
}
