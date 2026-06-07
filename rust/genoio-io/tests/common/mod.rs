use std::ops::Deref;
use std::path::Path;

use tempfile::TempDir;

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
