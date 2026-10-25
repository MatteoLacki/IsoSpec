use std::fmt;

/// What went wrong, by the C++ exception class the library threw (or, for
/// arguments the Rust side rejects before calling in, the class it would have
/// thrown).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum ErrorKind {
    /// Malformed input: a formula, sequence or CSV that does not parse, an
    /// unknown element or Unimod id, a probability outside its range, slices
    /// of mismatched lengths, an argument the algorithm cannot work with
    /// (`std::invalid_argument`).
    InvalidArgument,
    /// A number too large for the representation it must fit, such as an atom
    /// count past `i32` in a formula (`std::out_of_range`).
    OutOfRange,
    /// A size the library cannot allocate or address, such as a bin grid
    /// finer than doubles can resolve (`std::length_error`).
    Length,
    /// An operation whose result would not be finite, such as scaling an
    /// envelope's probabilities past the largest double (`std::overflow_error`).
    Overflow,
    /// Memory allocation failed (`std::bad_alloc`).
    OutOfMemory,
    /// Any other exception.
    Other,
}

/// An error from IsoSpec: its kind and the library's message.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Error {
    kind: ErrorKind,
    message: String,
}

impl Error {
    pub(crate) fn new(kind: ErrorKind, message: impl Into<String>) -> Self {
        Error {
            kind,
            message: message.into(),
        }
    }

    pub(crate) fn invalid(message: impl Into<String>) -> Self {
        Error::new(ErrorKind::InvalidArgument, message)
    }

    /// The kind of error.
    pub fn kind(&self) -> ErrorKind {
        self.kind
    }

    /// The message, as the library wrote it.
    pub fn message(&self) -> &str {
        &self.message
    }
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.message)
    }
}

impl std::error::Error for Error {}

// rust/shim/shim.h's trycatch prefixes every message with one letter naming
// the exception class.
impl From<cxx::Exception> for Error {
    fn from(e: cxx::Exception) -> Self {
        let what = e.what();
        let kind = match what.as_bytes().first() {
            Some(b'I') => ErrorKind::InvalidArgument,
            Some(b'R') => ErrorKind::OutOfRange,
            Some(b'L') => ErrorKind::Length,
            Some(b'O') => ErrorKind::Overflow,
            Some(b'M') => ErrorKind::OutOfMemory,
            _ => ErrorKind::Other,
        };
        Error::new(kind, what.get(1..).unwrap_or(""))
    }
}

/// `Result` with this crate's [`Error`].
pub type Result<T> = std::result::Result<T, Error>;
