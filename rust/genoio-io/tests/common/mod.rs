// pattern: Imperative Shell

use std::ops::Deref;
use std::path::Path;

use tempfile::TempDir;

#[expect(
    dead_code,
    reason = "shared integration-test helpers are compiled per test target"
)]
pub(crate) mod bcf;
pub(crate) mod bgen_output;
#[expect(
    dead_code,
    reason = "shared integration-test helpers are compiled per test target"
)]
pub(crate) mod dense;
pub(crate) mod output;
pub(crate) mod plink_output;
pub(crate) mod vcf_output;

pub struct TestDir {
    inner: TempDir,
}

impl Deref for TestDir {
    type Target = Path;

    fn deref(&self) -> &Self::Target {
        self.inner.path()
    }
}

pub fn unique_dir(name: &str) -> TestDir {
    let prefix = format!("genoio-{name}-");
    let inner = tempfile::Builder::new()
        .prefix(&prefix)
        .tempdir()
        .expect("test temp dir should be created");
    TestDir { inner }
}
