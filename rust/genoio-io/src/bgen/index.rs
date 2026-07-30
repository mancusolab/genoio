// pattern: Imperative Shell
//! Query same-path `.bgi` indexes for BGEN region pushdown.
//!
//! BGEN indexes are SQLite databases keyed by chromosome and position. This
//! module returns byte ranges in source order and validates that indexed reads
//! seek to the expected variant payloads.

use std::io::Seek;
use std::path::{Path, PathBuf};

use genoio_core::{GenoioError, RegionPredicate, VariantFilter};
use rusqlite::{params, Connection};

#[cfg(test)]
use super::session::BgenWorkProbe;
use crate::Result;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) struct BgenIndexRecord {
    pub(super) file_start_position: u64,
    size_in_bytes: u64,
}

pub(super) fn indexed_region_records(
    bgen: &Path,
    variant_filter: Option<&VariantFilter>,
) -> Result<Option<Vec<BgenIndexRecord>>> {
    #[cfg(test)]
    {
        indexed_region_records_inner(bgen, variant_filter, None)
    }
    #[cfg(not(test))]
    {
        indexed_region_records_inner(bgen, variant_filter)
    }
}

#[cfg(test)]
pub(super) fn indexed_region_records_with_probe(
    bgen: &Path,
    variant_filter: Option<&VariantFilter>,
    probe: &BgenWorkProbe,
) -> Result<Option<Vec<BgenIndexRecord>>> {
    indexed_region_records_inner(bgen, variant_filter, Some(probe))
}

fn indexed_region_records_inner(
    bgen: &Path,
    variant_filter: Option<&VariantFilter>,
    #[cfg(test)] probe: Option<&BgenWorkProbe>,
) -> Result<Option<Vec<BgenIndexRecord>>> {
    let Some(region) = variant_filter.and_then(VariantFilter::concrete_region_pushdown) else {
        return Ok(None);
    };
    let index_path = bgen_index_path(bgen);
    if !index_path.exists() {
        return Ok(None);
    }
    if !index_path.is_file() {
        return Err(GenoioError::invalid_source(
            &index_path,
            "bgen index path is not a file",
        ));
    }
    #[cfg(test)]
    {
        query_bgen_index_inner(&index_path, &region, probe).map(Some)
    }
    #[cfg(not(test))]
    {
        query_bgen_index(&index_path, &region).map(Some)
    }
}

fn bgen_index_path(bgen: &Path) -> PathBuf {
    let mut path = bgen.as_os_str().to_os_string();
    path.push(".bgi");
    PathBuf::from(path)
}

fn query_bgen_index(index_path: &Path, region: &RegionPredicate) -> Result<Vec<BgenIndexRecord>> {
    #[cfg(test)]
    {
        query_bgen_index_inner(index_path, region, None)
    }
    #[cfg(not(test))]
    {
        query_bgen_index_inner(index_path, region)
    }
}

fn query_bgen_index_inner(
    index_path: &Path,
    region: &RegionPredicate,
    #[cfg(test)] probe: Option<&BgenWorkProbe>,
) -> Result<Vec<BgenIndexRecord>> {
    let connection = Connection::open(index_path).map_err(|error| {
        GenoioError::invalid_source(index_path, format!("bgen index open error: {error}"))
    })?;
    #[cfg(test)]
    if let Some(probe) = probe {
        probe.record_index_open();
    }
    let mut statement = connection
        .prepare(
            // Preserve the BGEN source order contract even when records in the
            // indexed interval are not sorted by position.
            "SELECT file_start_position, size_in_bytes
             FROM Variant
             WHERE chromosome = ?1 AND position BETWEEN ?2 AND ?3
             ORDER BY file_start_position ASC, rowid ASC",
        )
        .map_err(|error| {
            GenoioError::invalid_source(
                index_path,
                format!("bgen index query prepare error: {error}"),
            )
        })?;
    let rows = statement
        .query_map(
            params![region.chrom, i64::from(region.start), i64::from(region.end)],
            |row| {
                let file_start_position = row.get::<_, i64>(0)?;
                let size_in_bytes = row.get::<_, i64>(1)?;
                Ok((file_start_position, size_in_bytes))
            },
        )
        .map_err(|error| {
            GenoioError::invalid_source(index_path, format!("bgen index query error: {error}"))
        })?;

    let mut records = Vec::new();
    for row in rows {
        let (file_start_position, size_in_bytes) = row.map_err(|error| {
            GenoioError::invalid_source(index_path, format!("bgen index row error: {error}"))
        })?;
        let file_start_position = u64::try_from(file_start_position).map_err(|_| {
            GenoioError::invalid_source(index_path, "bgen index file_start_position is negative")
        })?;
        let size_in_bytes = u64::try_from(size_in_bytes).map_err(|_| {
            GenoioError::invalid_source(index_path, "bgen index size_in_bytes is negative")
        })?;
        file_start_position
            .checked_add(size_in_bytes)
            .ok_or_else(|| {
                GenoioError::invalid_source(index_path, "bgen index byte range is out of range")
            })?;
        records.push(BgenIndexRecord {
            file_start_position,
            size_in_bytes,
        });
    }
    Ok(records)
}

pub(super) fn validate_index_record_consumed(
    reader: &mut impl Seek,
    bgen: &Path,
    index_record: &BgenIndexRecord,
) -> Result<()> {
    let consumed_end = reader.stream_position().map_err(|source| GenoioError::Io {
        path: bgen.to_path_buf(),
        source,
    })?;
    let expected_end = index_record
        .file_start_position
        .checked_add(index_record.size_in_bytes)
        .ok_or_else(|| {
            GenoioError::invalid_source(bgen, "bgen index byte range is out of range")
        })?;
    if consumed_end != expected_end {
        return Err(GenoioError::invalid_source(
            bgen,
            "bgen index byte range does not match decoded variant record",
        ));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use genoio_core::RegionPredicate;
    use rusqlite::params;

    use super::query_bgen_index;

    fn create_index_database() -> (tempfile::TempDir, std::path::PathBuf) {
        let directory = tempfile::tempdir().expect("test directory should be created");
        let path = directory.path().join("tiny.bgen.bgi");
        let connection = rusqlite::Connection::open(&path).expect("test index should open");
        connection
            .execute_batch(
                "CREATE TABLE Variant (
                    chromosome TEXT NOT NULL,
                    position INT NOT NULL,
                    file_start_position INT NOT NULL,
                    size_in_bytes INT NOT NULL
                );
                CREATE INDEX variant_source_offset
                    ON Variant(file_start_position ASC, size_in_bytes DESC);",
            )
            .expect("test index schema should be created");
        (directory, path)
    }

    fn region() -> RegionPredicate {
        RegionPredicate {
            chrom: "1".to_owned(),
            start: 1,
            end: 100,
        }
    }

    #[test]
    fn pbr_rust_bgen_003_tied_offsets_use_rowid_as_total_order() {
        let (_directory, path) = create_index_database();
        let connection = rusqlite::Connection::open(&path).expect("test index should open");
        for size in [11_i64, 22_i64] {
            connection
                .execute(
                    "INSERT INTO Variant (
                        chromosome, position, file_start_position, size_in_bytes
                     ) VALUES ('1', 10, 50, ?1)",
                    params![size],
                )
                .expect("test index row should insert");
        }
        drop(connection);

        let records = query_bgen_index(&path, &region())
            .expect("tied indexed offsets should query into owned records");

        assert_eq!(
            records
                .iter()
                .map(|record| record.size_in_bytes)
                .collect::<Vec<_>>(),
            vec![11, 22]
        );
    }

    #[test]
    fn pbr_rust_bgen_003_negative_index_scalars_are_rejected() {
        let (_directory, path) = create_index_database();
        let connection = rusqlite::Connection::open(&path).expect("test index should open");
        connection
            .execute(
                "INSERT INTO Variant (
                    chromosome, position, file_start_position, size_in_bytes
                 ) VALUES ('1', 10, -1, 20)",
                [],
            )
            .expect("test index row should insert");
        drop(connection);

        let error =
            query_bgen_index(&path, &region()).expect_err("negative indexed offset should fail");

        assert!(error
            .to_string()
            .contains("file_start_position is negative"));
    }
}
