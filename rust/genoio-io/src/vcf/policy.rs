// pattern: Mixed
// Reason: Text VCF routing combines pure filter-policy decisions with cheap
// filesystem checks for companion tabix/CSI indexes.

use std::ffi::OsString;
use std::path::{Path, PathBuf};

use genoio_core::{GenoioError, RegionPredicate, VariantFilter};

use crate::error::Result;

use super::is_compressed_vcf;

pub(super) fn read_text_vcf_with_optional_index<T, EmptyRead, IndexedRead, FullRead>(
    path: &Path,
    variant_filter: Option<&VariantFilter>,
    empty_read: EmptyRead,
    indexed_read: IndexedRead,
    full_read: FullRead,
) -> Result<T>
where
    EmptyRead: FnOnce() -> Result<T>,
    IndexedRead: FnOnce(&RegionPredicate) -> Result<T>,
    FullRead: FnOnce() -> Result<T>,
{
    if variant_filter.is_some_and(VariantFilter::is_always_false) {
        return empty_read();
    }

    // Concrete region predicates can use tabix/CSI chunks. Region expressions
    // that cannot be reduced to one interval stay on the permissive full scan.
    if let Some(region) = variant_filter.and_then(VariantFilter::concrete_region_pushdown) {
        if is_compressed_vcf(path) && has_vcf_index(path) {
            return indexed_read(&region);
        }
    }

    reject_unindexed_compressed_region(path, variant_filter)?;
    full_read()
}

pub(super) fn reject_unindexed_compressed_region(
    path: &Path,
    variant_filter: Option<&VariantFilter>,
) -> Result<()> {
    if variant_filter
        .and_then(VariantFilter::concrete_region_pushdown)
        .is_none()
        || !is_compressed_vcf(path)
    {
        return Ok(());
    }
    if has_vcf_index(path) {
        return Ok(());
    }
    Err(GenoioError::invalid_source(
        path,
        "region filter on compressed VCF requires an index",
    ))
}

pub(super) fn has_vcf_index(path: &Path) -> bool {
    companion_index_path(path, "tbi").exists() || companion_index_path(path, "csi").exists()
}

pub(super) fn companion_index_path(path: &Path, index_extension: &str) -> PathBuf {
    let mut raw = OsString::from(path);
    raw.push(".");
    raw.push(index_extension);
    PathBuf::from(raw)
}
