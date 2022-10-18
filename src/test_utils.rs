use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::{env, fmt, process};

// Yeast CHR01
pub const YEAST_CHR01_INDEX: &str = "test_data/yeast_chr01_index";

// Single contig
pub const TINY_INDEX: &str = "test_data/tiny_index";
pub const TINY_RC_INDEX: &str = "test_data/tiny-rc_index";

// Multi reference
pub const TINY_REFS_INDEX: &str = "test_data/tiny-multi-refs/tiny-multi-refs_index";

// Small txome
pub const SMALL_TXOME_DENSE_INDEX: &str = "test_data/small_txome_index";
pub const SMALL_TXOME_SPARSE_INDEX: &str = "test_data/small_txome_index_sparse";

pub fn assert_path_eq(p: PathBuf, s: &str) {
    assert_eq!(p.to_str().unwrap(), s);
}

pub fn to_abs_path<P: AsRef<Path>>(path: P) -> PathBuf {
    let workdir = env!("CARGO_MANIFEST_DIR");
    let workdir = Path::new(&workdir);
    workdir.join(path)
}

/// Temp file generateion taken from simple_sds
// Counter used for temporary file names.
static TEMP_FILE_COUNTER: AtomicUsize = AtomicUsize::new(0);

/// Returns a name for a temporary file using the provided name part.
///
/// # Examples
///
/// ```
/// use simple_sds::serialize;
///
/// let filename = serialize::temp_file_name("example");
/// assert!(filename.into_os_string().into_string().unwrap().contains("example"));
/// ```
pub fn temp_file_name(name_part: &str) -> PathBuf {
    let count = TEMP_FILE_COUNTER.fetch_add(1, Ordering::SeqCst);
    let mut buf = env::temp_dir();
    buf.push(format!("{}_{}_{}", name_part, process::id(), count));
    buf
}

/// Taken from itertools crate
/// https://docs.rs/itertools/0.5.9/src/itertools/lib.rs.html#1576-1599
///
/// Assert that two iterators produce equal sequences, with the same
/// semantics as *equal(a, b)*.
///
/// **Panics** on assertion failure with a message that shows the
/// two iteration elements.
///
/// ```ignore
/// assert_equal("exceed".split('c'), "excess".split('c'));
/// // ^PANIC: panicked at 'Failed assertion Some("eed") == Some("ess") for iteration 1',
/// ```
pub fn assert_iter_eq<I, J>(a: I, b: J)
where
    I: IntoIterator,
    J: IntoIterator,
    I::Item: fmt::Debug + PartialEq<J::Item>,
    J::Item: fmt::Debug,
{
    let mut ia = a.into_iter();
    let mut ib = b.into_iter();
    let mut i = 0;
    loop {
        match (ia.next(), ib.next()) {
            (None, None) => return,
            (Some(a), Some(b)) => {
                assert_eq!(a, b, "\nFailed for iteration {}\n", i,);
                i += 1;
            }
            (a, b) => {
                panic!(
                    "\nFailed for iteration {i}\n{a:#?}\n=/=\n{b:#?}\n",
                    i = i,
                    a = a,
                    b = b
                );
            }
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    #[test]
    #[should_panic]
    fn not_iter_eq() {
        assert_iter_eq([0, 1], [0, 0]);
    }

    #[test]
    fn iter_eq() {
        assert_iter_eq([0, 1], [0, 1]);
    }
}
