// pattern: Mixed (unavoidable)
// Reason: Format-local binary parsing is kept beside the filesystem entrypoint to match the
// existing reader module pattern in this crate.

use std::collections::HashSet;
use std::fs::{self, File};
use std::io::{Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};

use flate2::read::ZlibDecoder;
use genoio_core::{
    attach_variant_stats, compute_dosage_variant_stats, select_samples_source_order,
    transpose_variant_major_to_sample_major, DenseGenotypeMatrix, DenseSampleSelection,
    MetadataError, MetadataOutput, PartialFilterDecision, RegionPredicate, SampleRecord,
    SourceCapabilities, VariantFilter, VariantRecord, VariantWindow,
};
use rusqlite::{params, Connection};

use crate::Result;

const BGEN_MAGIC: &[u8; 4] = b"bgen";
const ZERO_MAGIC: &[u8; 4] = &[0, 0, 0, 0];
const MIN_HEADER_LENGTH: u32 = 20;
const SAMPLE_IDENTIFIER_FLAG: u32 = 1 << 31;
const DOSAGE_TOLERANCE: f32 = 1.0e-6;

pub fn read_bgen_metadata(bgen: &Path, sample: Option<&Path>) -> Result<MetadataOutput> {
    let mut reader = File::open(bgen).map_err(|source| MetadataError::Io {
        path: bgen.to_path_buf(),
        source,
    })?;
    let header = BgenHeader::read_from(&mut reader, bgen)?;
    header.validate(bgen)?;

    let samples = read_bgen_samples(&mut reader, bgen, sample, &header)?;

    reader
        .seek(SeekFrom::Start(u64::from(header.offset) + 4))
        .map_err(|source| MetadataError::Io {
            path: bgen.to_path_buf(),
            source,
        })?;
    let variants = read_layout2_variant_metadata(
        &mut reader,
        bgen,
        header.variant_count,
        header.sample_count,
        header.flags.compression,
    )?;

    Ok(MetadataOutput {
        samples,
        variants,
        capabilities: SourceCapabilities::genotype_only(),
    })
}

/// Read all retained BGEN biallelic diploid dosages as a dense matrix.
pub fn read_bgen_dosage_dense(
    bgen: &Path,
    sample: Option<&Path>,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
) -> Result<DenseGenotypeMatrix> {
    read_bgen_dosage_dense_windowed(bgen, sample, requested_samples, variant_filter, None, false)
}

/// Read retained BGEN biallelic diploid dosages as a dense matrix.
pub fn read_bgen_dosage_dense_windowed(
    bgen: &Path,
    sample: Option<&Path>,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    matrix_only: bool,
) -> Result<DenseGenotypeMatrix> {
    let mut reader = File::open(bgen).map_err(|source| MetadataError::Io {
        path: bgen.to_path_buf(),
        source,
    })?;
    let header = BgenHeader::read_from(&mut reader, bgen)?;
    header.validate(bgen)?;
    let all_samples = read_bgen_samples(&mut reader, bgen, sample, &header)?;
    let selection = select_samples_source_order(&all_samples, requested_samples, bgen)?;
    let mut diagnostics = selection.diagnostics.clone();
    if variant_filter.is_some_and(VariantFilter::is_always_false) {
        diagnostics.retained_variants = 0;
        return finish_bgen_dense_matrix(
            BgenDenseParts {
                n_samples: selection.samples.len(),
                n_variants: 0,
                values: Vec::new(),
                missing_mask: Vec::new(),
                samples: selection.samples,
                variants: Vec::new(),
                diagnostics,
            },
            matrix_only,
        );
    }
    if let Some(index_records) = indexed_region_records(bgen, variant_filter)? {
        let context = BgenIndexedReadContext {
            reader: &mut reader,
            bgen,
            header: &header,
            selection,
            diagnostics,
            variant_filter,
            variant_window,
            matrix_only,
        };
        return read_bgen_dosage_dense_indexed(context, &index_records);
    }

    reader
        .seek(SeekFrom::Start(u64::from(header.offset) + 4))
        .map_err(|source| MetadataError::Io {
            path: bgen.to_path_buf(),
            source,
        })?;

    let header_variant_count = usize::try_from(header.variant_count)
        .map_err(|_| MetadataError::parse(bgen, "bgen variant count is out of range"))?;
    let output_variant_capacity = variant_window.map_or(header_variant_count, |window| {
        window.len.min(header_variant_count)
    });
    let mut variants = Vec::with_capacity(output_variant_capacity);
    let mut variant_major_values =
        Vec::with_capacity(selection.samples.len() * output_variant_capacity);
    let mut variant_major_missing =
        Vec::with_capacity(selection.samples.len() * output_variant_capacity);
    let mut decode_buffers = DosageDecodeBuffers::default();
    let mut retained_index = 0_usize;
    let mut output_variant_count = 0_usize;

    for _ in 0..header.variant_count {
        let mut variant = read_layout2_variant_identifying_data(&mut reader, bgen)?;
        diagnostics.candidate_variants += 1;
        let partial_decision = variant_filter.map_or(PartialFilterDecision::Accept, |filter| {
            filter.partial_decision(&variant)
        });
        match partial_decision {
            PartialFilterDecision::Reject => {
                skip_layout2_probability_payload(&mut reader, bgen, header.flags.compression)?;
                diagnostics.dropped_metadata_variants += 1;
                continue;
            }
            PartialFilterDecision::Accept => {
                let include_in_window =
                    variant_window.is_none_or(|window| window.contains(retained_index));
                retained_index += 1;
                if !include_in_window {
                    if variant_window.is_some_and(|window| window.is_past(retained_index)) {
                        break;
                    }
                    skip_layout2_probability_payload(&mut reader, bgen, header.flags.compression)?;
                    continue;
                }
                read_layout2_probability_payload_into(
                    &mut reader,
                    bgen,
                    header.flags.compression,
                    &mut decode_buffers.payload,
                    &mut decode_buffers.compressed_payload,
                )?;
            }
            PartialFilterDecision::NeedGenotypes => {
                read_layout2_probability_payload_into(
                    &mut reader,
                    bgen,
                    header.flags.compression,
                    &mut decode_buffers.payload,
                    &mut decode_buffers.compressed_payload,
                )?;
                decode_buffered_dosage_values(
                    bgen,
                    header.sample_count,
                    &selection.source_indices,
                    &mut decode_buffers,
                )?;
                let stats = compute_dosage_variant_stats(
                    &decode_buffers.selected_values,
                    &decode_buffers.selected_missing,
                )?;
                if variant_filter.is_some_and(|filter| !filter.evaluate(&variant, Some(&stats))) {
                    diagnostics.dropped_genotype_variants += 1;
                    continue;
                }
                attach_variant_stats(&mut variant, stats);
                let include_in_window =
                    variant_window.is_none_or(|window| window.contains(retained_index));
                retained_index += 1;
                if !include_in_window {
                    if variant_window.is_some_and(|window| window.is_past(retained_index)) {
                        break;
                    }
                    continue;
                }
            }
        }

        if !matches!(partial_decision, PartialFilterDecision::NeedGenotypes) {
            decode_buffered_dosage_values(
                bgen,
                header.sample_count,
                &selection.source_indices,
                &mut decode_buffers,
            )?;
        }
        if !matrix_only {
            variants.push(variant);
        }
        variant_major_values.extend_from_slice(&decode_buffers.selected_values);
        variant_major_missing.extend_from_slice(&decode_buffers.selected_missing);
        output_variant_count += 1;
        if variant_window.is_some_and(|window| window.is_past(retained_index)) {
            break;
        }
    }

    let n_samples = selection.samples.len();
    let n_variants = output_variant_count;
    diagnostics.retained_variants = n_variants;
    let values =
        transpose_variant_major_to_sample_major(&variant_major_values, n_samples, n_variants);
    let missing_mask =
        transpose_variant_major_to_sample_major(&variant_major_missing, n_samples, n_variants);

    finish_bgen_dense_matrix(
        BgenDenseParts {
            n_samples,
            n_variants,
            values,
            missing_mask,
            samples: selection.samples,
            variants,
            diagnostics,
        },
        matrix_only,
    )
}

struct BgenDenseParts {
    n_samples: usize,
    n_variants: usize,
    values: Vec<f32>,
    missing_mask: Vec<bool>,
    samples: Vec<SampleRecord>,
    variants: Vec<VariantRecord>,
    diagnostics: genoio_core::DenseDiagnostics,
}

fn finish_bgen_dense_matrix(
    parts: BgenDenseParts,
    matrix_only: bool,
) -> Result<DenseGenotypeMatrix> {
    let BgenDenseParts {
        n_samples,
        n_variants,
        values,
        missing_mask,
        samples,
        variants,
        diagnostics,
    } = parts;
    if matrix_only {
        DenseGenotypeMatrix::new_matrix_only(
            n_samples,
            n_variants,
            values,
            missing_mask,
            diagnostics,
        )
    } else {
        DenseGenotypeMatrix::new(
            n_samples,
            n_variants,
            values,
            missing_mask,
            samples,
            variants,
            diagnostics,
        )
    }
}

struct BgenIndexedReadContext<'a> {
    reader: &'a mut File,
    bgen: &'a Path,
    header: &'a BgenHeader,
    selection: DenseSampleSelection,
    diagnostics: genoio_core::DenseDiagnostics,
    variant_filter: Option<&'a VariantFilter>,
    variant_window: Option<VariantWindow>,
    matrix_only: bool,
}

fn read_bgen_dosage_dense_indexed(
    context: BgenIndexedReadContext<'_>,
    index_records: &[BgenIndexRecord],
) -> Result<DenseGenotypeMatrix> {
    let BgenIndexedReadContext {
        reader,
        bgen,
        header,
        selection,
        mut diagnostics,
        variant_filter,
        variant_window,
        matrix_only,
    } = context;
    let output_variant_capacity = variant_window.map_or(index_records.len(), |window| {
        window.len.min(index_records.len())
    });
    let mut variants = Vec::with_capacity(output_variant_capacity);
    let mut variant_major_values =
        Vec::with_capacity(selection.samples.len() * output_variant_capacity);
    let mut variant_major_missing =
        Vec::with_capacity(selection.samples.len() * output_variant_capacity);
    let mut decode_buffers = DosageDecodeBuffers::default();
    let mut retained_index = 0_usize;
    let mut output_variant_count = 0_usize;

    for index_record in index_records {
        if variant_window.is_some_and(|window| window.is_past(retained_index)) {
            break;
        }
        reader
            .seek(SeekFrom::Start(index_record.file_start_position))
            .map_err(|source| MetadataError::Io {
                path: bgen.to_path_buf(),
                source,
            })?;
        let mut variant = read_layout2_variant_identifying_data(reader, bgen)?;
        diagnostics.candidate_variants += 1;
        let partial_decision = variant_filter.map_or(PartialFilterDecision::Accept, |filter| {
            filter.partial_decision(&variant)
        });
        match partial_decision {
            PartialFilterDecision::Reject => {
                skip_layout2_probability_payload(reader, bgen, header.flags.compression)?;
                validate_index_record_consumed(reader, bgen, index_record)?;
                diagnostics.dropped_metadata_variants += 1;
                continue;
            }
            PartialFilterDecision::Accept => {
                let include_in_window =
                    variant_window.is_none_or(|window| window.contains(retained_index));
                retained_index += 1;
                if !include_in_window {
                    skip_layout2_probability_payload(reader, bgen, header.flags.compression)?;
                    validate_index_record_consumed(reader, bgen, index_record)?;
                    continue;
                }
                read_layout2_probability_payload_into(
                    reader,
                    bgen,
                    header.flags.compression,
                    &mut decode_buffers.payload,
                    &mut decode_buffers.compressed_payload,
                )?;
            }
            PartialFilterDecision::NeedGenotypes => {
                read_layout2_probability_payload_into(
                    reader,
                    bgen,
                    header.flags.compression,
                    &mut decode_buffers.payload,
                    &mut decode_buffers.compressed_payload,
                )?;
                decode_buffered_dosage_values(
                    bgen,
                    header.sample_count,
                    &selection.source_indices,
                    &mut decode_buffers,
                )?;
                let stats = compute_dosage_variant_stats(
                    &decode_buffers.selected_values,
                    &decode_buffers.selected_missing,
                )?;
                if variant_filter.is_some_and(|filter| !filter.evaluate(&variant, Some(&stats))) {
                    validate_index_record_consumed(reader, bgen, index_record)?;
                    diagnostics.dropped_genotype_variants += 1;
                    continue;
                }
                attach_variant_stats(&mut variant, stats);
                let include_in_window =
                    variant_window.is_none_or(|window| window.contains(retained_index));
                retained_index += 1;
                if !include_in_window {
                    validate_index_record_consumed(reader, bgen, index_record)?;
                    continue;
                }
            }
        }

        if !matches!(partial_decision, PartialFilterDecision::NeedGenotypes) {
            decode_buffered_dosage_values(
                bgen,
                header.sample_count,
                &selection.source_indices,
                &mut decode_buffers,
            )?;
        }
        validate_index_record_consumed(reader, bgen, index_record)?;
        if !matrix_only {
            variants.push(variant);
        }
        variant_major_values.extend_from_slice(&decode_buffers.selected_values);
        variant_major_missing.extend_from_slice(&decode_buffers.selected_missing);
        output_variant_count += 1;
    }

    let n_samples = selection.samples.len();
    let n_variants = output_variant_count;
    diagnostics.retained_variants = n_variants;
    let values =
        transpose_variant_major_to_sample_major(&variant_major_values, n_samples, n_variants);
    let missing_mask =
        transpose_variant_major_to_sample_major(&variant_major_missing, n_samples, n_variants);

    finish_bgen_dense_matrix(
        BgenDenseParts {
            n_samples,
            n_variants,
            values,
            missing_mask,
            samples: selection.samples,
            variants,
            diagnostics,
        },
        matrix_only,
    )
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct BgenIndexRecord {
    file_start_position: u64,
    size_in_bytes: u64,
}

fn indexed_region_records(
    bgen: &Path,
    variant_filter: Option<&VariantFilter>,
) -> Result<Option<Vec<BgenIndexRecord>>> {
    let Some(region) = variant_filter.and_then(VariantFilter::concrete_region_pushdown) else {
        return Ok(None);
    };
    let index_path = bgen_index_path(bgen);
    if !index_path.exists() {
        return Ok(None);
    }
    if !index_path.is_file() {
        return Err(MetadataError::parse(
            &index_path,
            "bgen index path is not a file",
        ));
    }
    query_bgen_index(&index_path, &region).map(Some)
}

fn bgen_index_path(bgen: &Path) -> PathBuf {
    let mut path = bgen.as_os_str().to_os_string();
    path.push(".bgi");
    PathBuf::from(path)
}

fn query_bgen_index(index_path: &Path, region: &RegionPredicate) -> Result<Vec<BgenIndexRecord>> {
    let connection = Connection::open(index_path).map_err(|error| {
        MetadataError::parse(index_path, format!("bgen index open error: {error}"))
    })?;
    let mut statement = connection
        .prepare(
            // Preserve the BGEN source order contract even when records in the
            // indexed interval are not sorted by position.
            "SELECT file_start_position, size_in_bytes
             FROM Variant
             WHERE chromosome = ?1 AND position BETWEEN ?2 AND ?3
             ORDER BY file_start_position",
        )
        .map_err(|error| {
            MetadataError::parse(
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
            MetadataError::parse(index_path, format!("bgen index query error: {error}"))
        })?;

    let mut records = Vec::new();
    for row in rows {
        let (file_start_position, size_in_bytes) = row.map_err(|error| {
            MetadataError::parse(index_path, format!("bgen index row error: {error}"))
        })?;
        records.push(BgenIndexRecord {
            file_start_position: u64::try_from(file_start_position).map_err(|_| {
                MetadataError::parse(index_path, "bgen index file_start_position is negative")
            })?,
            size_in_bytes: u64::try_from(size_in_bytes).map_err(|_| {
                MetadataError::parse(index_path, "bgen index size_in_bytes is negative")
            })?,
        });
    }
    Ok(records)
}

fn validate_index_record_consumed(
    reader: &mut File,
    bgen: &Path,
    index_record: &BgenIndexRecord,
) -> Result<()> {
    let consumed_end = reader
        .stream_position()
        .map_err(|source| MetadataError::Io {
            path: bgen.to_path_buf(),
            source,
        })?;
    let expected_end = index_record
        .file_start_position
        .checked_add(index_record.size_in_bytes)
        .ok_or_else(|| MetadataError::parse(bgen, "bgen index byte range is out of range"))?;
    if consumed_end != expected_end {
        return Err(MetadataError::parse(
            bgen,
            "bgen index byte range does not match decoded variant record",
        ));
    }
    Ok(())
}

#[derive(Default)]
struct DosageDecodeBuffers {
    payload: Vec<u8>,
    compressed_payload: Vec<u8>,
    selected_values: Vec<f32>,
    selected_missing: Vec<bool>,
}

fn decode_buffered_dosage_values(
    bgen: &Path,
    sample_count: u32,
    source_indices: &[usize],
    buffers: &mut DosageDecodeBuffers,
) -> Result<()> {
    let DosageDecodeBuffers {
        payload,
        selected_values,
        selected_missing,
        ..
    } = buffers;
    let decoded = DecodedDosageVariant::decode(bgen, payload, sample_count, 2)?;
    decoded.debug_assert_supported_subset();
    decoded.decode_selected_source_order(
        bgen,
        source_indices,
        selected_values,
        selected_missing,
    )?;
    Ok(())
}

fn read_bgen_samples(
    reader: &mut impl Read,
    bgen: &Path,
    sample: Option<&Path>,
    header: &BgenHeader,
) -> Result<Vec<SampleRecord>> {
    if header.flags.has_sample_identifiers {
        read_sample_identifier_block(reader, bgen, header.sample_count)
    } else if let Some(sample) = sample {
        read_companion_sample_file(sample, header.sample_count)
    } else {
        Err(MetadataError::parse(
            bgen,
            "bgen sample identifiers require embedded identifiers or a companion sample path",
        ))
    }
}

struct BgenHeader {
    offset: u32,
    header_length: u32,
    variant_count: u32,
    sample_count: u32,
    flags: BgenFlags,
}

impl BgenHeader {
    fn read_from(reader: &mut impl Read, path: &Path) -> Result<Self> {
        let offset = read_u32_le(reader, path)?;
        let header_length = read_u32_le(reader, path)?;
        let variant_count = read_u32_le(reader, path)?;
        let sample_count = read_u32_le(reader, path)?;

        let mut magic = [0_u8; 4];
        reader
            .read_exact(&mut magic)
            .map_err(|source| MetadataError::Io {
                path: path.to_path_buf(),
                source,
            })?;
        if &magic != BGEN_MAGIC && &magic != ZERO_MAGIC {
            return Err(MetadataError::parse(path, "invalid bgen magic bytes"));
        }

        let free_data_length = header_length
            .checked_sub(MIN_HEADER_LENGTH)
            .ok_or_else(|| {
                MetadataError::parse(path, "bgen header length is smaller than 20 bytes")
            })?;
        skip_exact(reader, path, u64::from(free_data_length))?;
        let flags = BgenFlags::from_raw(read_u32_le(reader, path)?);
        Ok(Self {
            offset,
            header_length,
            variant_count,
            sample_count,
            flags,
        })
    }

    fn validate(&self, path: &Path) -> Result<()> {
        if self.header_length < MIN_HEADER_LENGTH {
            return Err(MetadataError::parse(
                path,
                "bgen header length is smaller than 20 bytes",
            ));
        }
        if self.header_length > self.offset {
            return Err(MetadataError::parse(
                path,
                "bgen header length exceeds variant data offset",
            ));
        }
        if self.flags.layout != BgenLayout::Layout2 {
            return Err(MetadataError::parse(
                path,
                "bgen metadata parsing requires layout 2",
            ));
        }
        if self.flags.compression == BgenCompression::Reserved {
            return Err(MetadataError::parse(
                path,
                "bgen compression value is reserved",
            ));
        }
        Ok(())
    }
}

struct BgenFlags {
    compression: BgenCompression,
    layout: BgenLayout,
    has_sample_identifiers: bool,
}

impl BgenFlags {
    fn from_raw(raw: u32) -> Self {
        Self {
            compression: BgenCompression::from_raw(raw & 0b11),
            layout: BgenLayout::from_raw((raw >> 2) & 0b1111),
            has_sample_identifiers: raw & SAMPLE_IDENTIFIER_FLAG != 0,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum BgenCompression {
    None,
    Zlib,
    Zstd,
    Reserved,
}

impl BgenCompression {
    fn from_raw(raw: u32) -> Self {
        match raw {
            0 => Self::None,
            1 => Self::Zlib,
            2 => Self::Zstd,
            _ => Self::Reserved,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum BgenLayout {
    Layout1,
    Layout2,
    Reserved,
}

impl BgenLayout {
    fn from_raw(raw: u32) -> Self {
        match raw {
            1 => Self::Layout1,
            2 => Self::Layout2,
            _ => Self::Reserved,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum BgenPhase {
    Unphased,
    Phased,
}

impl BgenPhase {
    fn from_raw(path: &Path, raw: u8) -> Result<Self> {
        match raw {
            0 => Ok(Self::Unphased),
            1 => Ok(Self::Phased),
            _ => Err(MetadataError::parse(
                path,
                "unsupported bgen phased probability value; expected 0 or 1",
            )),
        }
    }
}

fn read_sample_identifier_block(
    reader: &mut impl Read,
    path: &Path,
    expected_sample_count: u32,
) -> Result<Vec<SampleRecord>> {
    let block_length = read_u32_le(reader, path)?;
    if block_length < 8 {
        return Err(MetadataError::parse(
            path,
            "bgen sample identifiers block is too short",
        ));
    }

    let sample_count = read_u32_le(reader, path)?;
    if sample_count != expected_sample_count {
        return Err(MetadataError::parse(
            path,
            "bgen sample identifiers count does not match header sample count",
        ));
    }
    let body_len = usize::try_from(block_length - 8).map_err(|_| {
        MetadataError::parse(path, "bgen sample identifiers block length is out of range")
    })?;
    let body = read_exact_vec(reader, path, body_len)?;
    let mut body_reader = body.as_slice();

    let mut records = Vec::with_capacity(usize::try_from(sample_count).map_err(|_| {
        MetadataError::parse(path, "bgen sample identifiers count is out of range")
    })?);
    let mut seen = HashSet::with_capacity(records.capacity());
    for _ in 0..sample_count {
        let sample_id = read_len_prefixed_string_u16(&mut body_reader, path, "sample identifier")?;
        if sample_id.is_empty() {
            return Err(MetadataError::parse(
                path,
                "bgen sample identifier is empty",
            ));
        }
        if !seen.insert(sample_id.clone()) {
            return Err(MetadataError::parse(
                path,
                format!("bgen duplicate sample identifier: {sample_id}"),
            ));
        }
        records.push(sample_record(sample_id));
    }

    if !body_reader.is_empty() {
        return Err(MetadataError::parse(
            path,
            "bgen sample identifiers block length does not match decoded sample identifiers",
        ));
    }

    Ok(records)
}

fn read_companion_sample_file(
    path: &Path,
    expected_sample_count: u32,
) -> Result<Vec<SampleRecord>> {
    let contents = fs::read_to_string(path).map_err(|source| MetadataError::Io {
        path: path.to_path_buf(),
        source,
    })?;
    let capacity = usize::try_from(expected_sample_count)
        .map_err(|_| MetadataError::parse(path, "bgen sample count is out of range"))?;
    let mut records = Vec::with_capacity(capacity);
    let mut seen = HashSet::with_capacity(capacity);

    for line in contents.lines().skip(2) {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let sample_id = line.split_whitespace().next().ok_or_else(|| {
            MetadataError::parse(path, "bgen companion sample identifier is empty")
        })?;
        if sample_id.is_empty() {
            return Err(MetadataError::parse(
                path,
                "bgen companion sample identifier is empty",
            ));
        }
        if !seen.insert(sample_id.to_owned()) {
            return Err(MetadataError::parse(
                path,
                format!("bgen duplicate sample identifier: {sample_id}"),
            ));
        }
        records.push(sample_record(sample_id.to_owned()));
    }

    if records.len() != capacity {
        return Err(MetadataError::parse(
            path,
            format!(
                "bgen companion sample count does not match header sample count: expected {capacity}, found {}",
                records.len()
            ),
        ));
    }

    Ok(records)
}

fn sample_record(iid: String) -> SampleRecord {
    SampleRecord {
        fid: None,
        iid,
        father: None,
        mother: None,
        sex: None,
        phenotype: None,
        source_sample_index: None,
        haplotype_index: None,
    }
}

fn read_layout2_variant_metadata(
    reader: &mut impl Read,
    path: &Path,
    variant_count: u32,
    sample_count: u32,
    compression: BgenCompression,
) -> Result<Vec<VariantRecord>> {
    let mut variants = Vec::with_capacity(
        usize::try_from(variant_count)
            .map_err(|_| MetadataError::parse(path, "bgen variant count is out of range"))?,
    );

    for _ in 0..variant_count {
        let variant = read_layout2_variant_identifying_data(reader, path)?;
        skip_layout2_probability_block(reader, path, sample_count, 2, compression)?;
        variants.push(variant);
    }

    if variants.len() != usize::try_from(variant_count).unwrap_or(usize::MAX) {
        return Err(MetadataError::parse(
            path,
            "bgen parsed variant count does not match header variant count",
        ));
    }

    Ok(variants)
}

fn read_layout2_variant_identifying_data(
    reader: &mut impl Read,
    path: &Path,
) -> Result<VariantRecord> {
    let id = read_len_prefixed_string_u16(reader, path, "variant id")?;
    let rsid = read_len_prefixed_string_u16(reader, path, "variant rsid")?;
    let chrom = read_len_prefixed_string_u16(reader, path, "variant chromosome")?;
    let pos = read_u32_le(reader, path)?;
    let allele_count = read_u16_le(reader, path)?;
    if allele_count != 2 {
        return Err(MetadataError::parse(
            path,
            "unsupported bgen multiallelic variant metadata; only biallelic records are supported",
        ));
    }

    let mut alleles = Vec::with_capacity(usize::from(allele_count));
    for _ in 0..allele_count {
        alleles.push(read_len_prefixed_string_u32(
            reader,
            path,
            "variant allele",
        )?);
    }
    let a0 = alleles[0].clone();
    let a1 = alleles[1].clone();
    let id = if rsid.is_empty() { id } else { rsid };

    Ok(VariantRecord {
        chrom,
        pos,
        id,
        a0: a0.clone(),
        a1: a1.clone(),
        ref_allele: Some(a0.clone()),
        alt_allele: Some(a1.clone()),
        source_a0: a0,
        source_a1: a1,
        flipped: false,
        qual: None,
        af: None,
        maf: None,
        mac: None,
        missing_rate: None,
        n_called: None,
    })
}

fn skip_layout2_probability_block(
    reader: &mut impl Read,
    path: &Path,
    expected_sample_count: u32,
    variant_allele_count: u16,
    compression: BgenCompression,
) -> Result<()> {
    let payload = read_layout2_probability_payload(reader, path, compression)?;
    let decoded =
        DecodedDosageVariant::decode(path, &payload, expected_sample_count, variant_allele_count)?;
    decoded.debug_assert_supported_subset();
    let mut values = Vec::new();
    let mut missing = Vec::new();
    decoded.decode_source_order(path, &mut values, &mut missing)?;
    Ok(())
}

fn read_layout2_probability_payload(
    reader: &mut impl Read,
    path: &Path,
    compression: BgenCompression,
) -> Result<Vec<u8>> {
    let mut payload = Vec::new();
    let mut compressed_payload = Vec::new();
    read_layout2_probability_payload_into(
        reader,
        path,
        compression,
        &mut payload,
        &mut compressed_payload,
    )?;
    Ok(payload)
}

fn read_layout2_probability_payload_into(
    reader: &mut impl Read,
    path: &Path,
    compression: BgenCompression,
    payload: &mut Vec<u8>,
    compressed_payload: &mut Vec<u8>,
) -> Result<()> {
    let block_length = read_u32_le(reader, path)?;
    match compression {
        BgenCompression::None => {
            let payload_length = usize::try_from(block_length).map_err(|_| {
                MetadataError::parse(path, "bgen uncompressed probability block is out of range")
            })?;
            payload.clear();
            payload.resize(payload_length, 0);
            reader
                .read_exact(payload)
                .map_err(|source| MetadataError::Io {
                    path: path.to_path_buf(),
                    source,
                })?;
            Ok(())
        }
        BgenCompression::Zlib | BgenCompression::Zstd => {
            if block_length < 4 {
                return Err(MetadataError::parse(
                    path,
                    "bgen compressed probability block length is smaller than decompressed length prefix",
                ));
            }
            let decompressed_block_length = read_u32_le(reader, path)?;
            let compressed_payload_length = usize::try_from(block_length - 4).map_err(|_| {
                MetadataError::parse(path, "bgen compressed probability block is out of range")
            })?;
            compressed_payload.clear();
            compressed_payload.resize(compressed_payload_length, 0);
            reader
                .read_exact(compressed_payload)
                .map_err(|source| MetadataError::Io {
                    path: path.to_path_buf(),
                    source,
                })?;
            decompress_probability_block_into(
                path,
                compression,
                compressed_payload,
                decompressed_block_length,
                payload,
            )
        }
        BgenCompression::Reserved => Err(MetadataError::parse(
            path,
            "bgen compression value is reserved",
        )),
    }
}

fn skip_layout2_probability_payload(
    reader: &mut impl Read,
    path: &Path,
    compression: BgenCompression,
) -> Result<()> {
    let block_length = read_u32_le(reader, path)?;
    match compression {
        BgenCompression::None | BgenCompression::Zlib | BgenCompression::Zstd => {
            skip_exact(reader, path, u64::from(block_length))
        }
        BgenCompression::Reserved => Err(MetadataError::parse(
            path,
            "bgen compression value is reserved",
        )),
    }
}

#[derive(Debug, Clone)]
struct Layout2ProbabilityHeader {
    sample_count: u32,
    allele_count: u16,
    min_ploidy: u8,
    max_ploidy: u8,
    sample_ploidies: Vec<u8>,
    non_missing_sample_count: u32,
    phase: BgenPhase,
    bit_depth: u8,
    byte_len: usize,
}

impl Layout2ProbabilityHeader {
    fn decode(
        path: &Path,
        payload: &[u8],
        expected_sample_count: u32,
        variant_allele_count: u16,
    ) -> Result<Self> {
        let fixed_header_length = Self::fixed_header_length(path, expected_sample_count)?;
        if payload.len() < fixed_header_length {
            return Err(MetadataError::parse(
                path,
                "bgen uncompressed probability block is shorter than the layout 2 header",
            ));
        }

        let reader = &mut &payload[..fixed_header_length];
        let sample_count = read_u32_le(reader, path)?;
        if sample_count != expected_sample_count {
            return Err(MetadataError::parse(
                path,
                "bgen probability block sample count does not match header sample count",
            ));
        }

        let allele_count = read_u16_le(reader, path)?;
        if allele_count != variant_allele_count {
            return Err(MetadataError::parse(
                path,
                "bgen probability block allele count does not match variant allele count",
            ));
        }
        if allele_count != 2 {
            return Err(MetadataError::parse(
                path,
                "unsupported bgen multiallelic probability block; only biallelic records are supported",
            ));
        }

        let min_ploidy = read_u8(reader, path)?;
        let max_ploidy = read_u8(reader, path)?;
        if min_ploidy != 2 || max_ploidy != 2 {
            return Err(MetadataError::parse(
                path,
                "unsupported bgen variable-ploidy probability block; only diploid records are supported",
            ));
        }

        let sample_count_usize = usize::try_from(expected_sample_count)
            .map_err(|_| MetadataError::parse(path, "bgen sample count is out of range"))?;
        let mut sample_ploidies = Vec::with_capacity(sample_count_usize);
        let mut non_missing_sample_count = 0_u32;
        for _ in 0..expected_sample_count {
            let ploidy_byte = read_u8(reader, path)?;
            let is_missing = ploidy_byte & 0b1000_0000 != 0;
            let ploidy = ploidy_byte & 0b0011_1111;
            if !is_missing {
                if ploidy != 2 {
                    return Err(MetadataError::parse(
                        path,
                        "unsupported bgen variable-ploidy probability block; only diploid records are supported",
                    ));
                }
                non_missing_sample_count =
                    non_missing_sample_count.checked_add(1).ok_or_else(|| {
                        MetadataError::parse(path, "bgen non-missing sample count is out of range")
                    })?;
            }
            sample_ploidies.push(ploidy_byte);
        }

        let phase = BgenPhase::from_raw(path, read_u8(reader, path)?)?;

        let bit_depth = read_u8(reader, path)?;
        if !(1..=32).contains(&bit_depth) {
            return Err(MetadataError::parse(
                path,
                "unsupported bgen probability bit depth; expected 1..=32",
            ));
        }

        Ok(Self {
            sample_count,
            allele_count,
            min_ploidy,
            max_ploidy,
            sample_ploidies,
            non_missing_sample_count,
            phase,
            bit_depth,
            byte_len: fixed_header_length,
        })
    }

    fn fixed_header_length(path: &Path, sample_count: u32) -> Result<usize> {
        let sample_count = usize::try_from(sample_count)
            .map_err(|_| MetadataError::parse(path, "bgen sample count is out of range"))?;
        10_usize
            .checked_add(sample_count)
            .ok_or_else(|| MetadataError::parse(path, "bgen probability header is out of range"))
    }

    fn required_packed_probability_bytes(&self, path: &Path) -> Result<usize> {
        let stored_probability_count = match self.phase {
            // For the supported biallelic diploid subset, unphased records
            // store P(A0/A0) and P(A0/A1); P(A1/A1) is inferred.
            BgenPhase::Unphased => 2,
            // Phased records store one non-last allele probability per
            // haplotype, so diploid biallelic records also store two values.
            BgenPhase::Phased => 2,
        };
        let bits = u64::from(self.non_missing_sample_count)
            .checked_mul(stored_probability_count)
            .and_then(|value| value.checked_mul(u64::from(self.bit_depth)))
            .ok_or_else(|| {
                MetadataError::parse(path, "bgen packed probability bit count is out of range")
            })?;
        let bytes = bits.div_ceil(8);
        usize::try_from(bytes).map_err(|_| {
            MetadataError::parse(path, "bgen packed probability bytes are out of range")
        })
    }
}

#[derive(Debug, Clone)]
struct DecodedDosageVariant<'a> {
    header: Layout2ProbabilityHeader,
    packed_probabilities: &'a [u8],
}

impl<'a> DecodedDosageVariant<'a> {
    fn decode(
        path: &Path,
        payload: &'a [u8],
        expected_sample_count: u32,
        variant_allele_count: u16,
    ) -> Result<Self> {
        let header = Layout2ProbabilityHeader::decode(
            path,
            payload,
            expected_sample_count,
            variant_allele_count,
        )?;
        let packed_probabilities = &payload[header.byte_len..];
        let required_packed_len = header.required_packed_probability_bytes(path)?;
        if packed_probabilities.len() < required_packed_len {
            return Err(MetadataError::parse(
                path,
                "bgen probability block is truncated; packed probabilities are shorter than declared non-missing samples",
            ));
        }
        if packed_probabilities.len() > required_packed_len {
            return Err(MetadataError::parse(
                path,
                "bgen probability block has trailing packed probability bytes",
            ));
        }

        Ok(Self {
            header,
            packed_probabilities,
        })
    }

    fn debug_assert_supported_subset(&self) {
        debug_assert_eq!(
            self.header.sample_count as usize,
            self.header.sample_ploidies.len()
        );
        debug_assert_eq!(self.header.allele_count, 2);
        debug_assert_eq!(self.header.min_ploidy, 2);
        debug_assert_eq!(self.header.max_ploidy, 2);
        debug_assert!((1..=32).contains(&self.header.bit_depth));
        debug_assert!(
            !self.packed_probabilities.is_empty() || self.header.non_missing_sample_count == 0
        );
    }

    fn decode_source_order(
        &self,
        path: &Path,
        values: &mut Vec<f32>,
        missing: &mut Vec<bool>,
    ) -> Result<()> {
        let sample_count = usize::try_from(self.header.sample_count)
            .map_err(|_| MetadataError::parse(path, "bgen sample count is out of range"))?;
        values.clear();
        missing.clear();
        values.reserve(sample_count);
        missing.reserve(sample_count);

        let mut bit_reader = LittleEndianBitReader::new(self.packed_probabilities);
        for &ploidy_byte in &self.header.sample_ploidies {
            let is_missing = ploidy_byte & 0b1000_0000 != 0;
            if is_missing {
                values.push(0.0);
                missing.push(true);
                continue;
            }

            let first_raw = bit_reader.read_u32(path, self.header.bit_depth)?;
            let second_raw = bit_reader.read_u32(path, self.header.bit_depth)?;
            let value = match self.header.phase {
                BgenPhase::Unphased => {
                    decode_unphased_a1_dosage(path, self.header.bit_depth, first_raw, second_raw)?
                }
                BgenPhase::Phased => {
                    decode_phased_a1_dosage(self.header.bit_depth, first_raw, second_raw)
                }
            };
            values.push(value);
            missing.push(false);
        }

        Ok(())
    }

    fn decode_selected_source_order(
        &self,
        path: &Path,
        source_indices: &[usize],
        values: &mut Vec<f32>,
        missing: &mut Vec<bool>,
    ) -> Result<()> {
        values.clear();
        missing.clear();
        values.reserve(source_indices.len());
        missing.reserve(source_indices.len());

        let mut selected_cursor = 0_usize;
        let mut bit_reader = LittleEndianBitReader::new(self.packed_probabilities);
        for (sample_index, &ploidy_byte) in self.header.sample_ploidies.iter().enumerate() {
            let is_selected = source_indices
                .get(selected_cursor)
                .is_some_and(|&source_index| source_index == sample_index);
            let is_missing = ploidy_byte & 0b1000_0000 != 0;
            if is_missing {
                if is_selected {
                    values.push(0.0);
                    missing.push(true);
                    selected_cursor += 1;
                }
                continue;
            }

            let first_raw = bit_reader.read_u32(path, self.header.bit_depth)?;
            let second_raw = bit_reader.read_u32(path, self.header.bit_depth)?;
            if is_selected {
                let value = match self.header.phase {
                    BgenPhase::Unphased => decode_unphased_a1_dosage(
                        path,
                        self.header.bit_depth,
                        first_raw,
                        second_raw,
                    )?,
                    BgenPhase::Phased => {
                        decode_phased_a1_dosage(self.header.bit_depth, first_raw, second_raw)
                    }
                };
                values.push(value);
                missing.push(false);
                selected_cursor += 1;
            }
        }
        debug_assert_eq!(selected_cursor, source_indices.len());

        Ok(())
    }
}

struct LittleEndianBitReader<'a> {
    bytes: &'a [u8],
    bit_offset: usize,
}

impl<'a> LittleEndianBitReader<'a> {
    fn new(bytes: &'a [u8]) -> Self {
        Self {
            bytes,
            bit_offset: 0,
        }
    }

    fn read_u32(&mut self, path: &Path, bit_count: u8) -> Result<u32> {
        debug_assert!((1..=32).contains(&bit_count));
        let available_bits = self.bytes.len().checked_mul(8).ok_or_else(|| {
            MetadataError::parse(path, "bgen packed probability bit count is out of range")
        })?;
        let end_bit = self
            .bit_offset
            .checked_add(usize::from(bit_count))
            .ok_or_else(|| {
                MetadataError::parse(path, "bgen packed probability bit count is out of range")
            })?;
        if end_bit > available_bits {
            return Err(MetadataError::parse(
                path,
                "bgen packed probability bits are truncated",
            ));
        }

        if self.bit_offset.is_multiple_of(8) {
            let byte_offset = self.bit_offset / 8;
            let value = match bit_count {
                8 => u32::from(self.bytes[byte_offset]),
                16 => u32::from(u16::from_le_bytes([
                    self.bytes[byte_offset],
                    self.bytes[byte_offset + 1],
                ])),
                32 => u32::from_le_bytes([
                    self.bytes[byte_offset],
                    self.bytes[byte_offset + 1],
                    self.bytes[byte_offset + 2],
                    self.bytes[byte_offset + 3],
                ]),
                _ => {
                    let mut value = 0_u32;
                    for output_bit in 0..bit_count {
                        let source_bit = self.bit_offset + usize::from(output_bit);
                        let byte = self.bytes[source_bit / 8];
                        let bit = (byte >> (source_bit % 8)) & 1;
                        value |= u32::from(bit) << output_bit;
                    }
                    self.bit_offset = end_bit;
                    return Ok(value);
                }
            };
            self.bit_offset = end_bit;
            return Ok(value);
        }

        let mut value = 0_u32;
        for output_bit in 0..bit_count {
            let source_bit = self.bit_offset + usize::from(output_bit);
            let byte = self.bytes[source_bit / 8];
            let bit = (byte >> (source_bit % 8)) & 1;
            value |= u32::from(bit) << output_bit;
        }
        self.bit_offset = end_bit;
        Ok(value)
    }
}

fn decode_unphased_a1_dosage(
    path: &Path,
    bit_depth: u8,
    p_aa_raw: u32,
    p_ab_raw: u32,
) -> Result<f32> {
    let denominator = probability_denominator(bit_depth);
    let p_aa = p_aa_raw as f32 / denominator;
    let p_ab = p_ab_raw as f32 / denominator;
    let p_bb = 1.0 - p_aa - p_ab;
    let dosage = p_ab + 2.0 * p_bb;

    if p_bb < -DOSAGE_TOLERANCE || !(-DOSAGE_TOLERANCE..=2.0 + DOSAGE_TOLERANCE).contains(&dosage) {
        return Err(MetadataError::parse(
            path,
            "bgen malformed probability values produce invalid a1 dosage",
        ));
    }

    Ok(dosage.clamp(0.0, 2.0))
}

fn decode_phased_a1_dosage(bit_depth: u8, p_hap0_a0_raw: u32, p_hap1_a0_raw: u32) -> f32 {
    let denominator = probability_denominator(bit_depth);
    let p_hap0_a0 = p_hap0_a0_raw as f32 / denominator;
    let p_hap1_a0 = p_hap1_a0_raw as f32 / denominator;
    (2.0 - p_hap0_a0 - p_hap1_a0).clamp(0.0, 2.0)
}

fn probability_denominator(bit_depth: u8) -> f32 {
    ((1_u64 << bit_depth) - 1) as f32
}

fn decompress_probability_block_into(
    path: &Path,
    compression: BgenCompression,
    compressed_payload: &[u8],
    expected_decompressed_len: u32,
    decompressed: &mut Vec<u8>,
) -> Result<()> {
    let capacity = usize::try_from(expected_decompressed_len).map_err(|_| {
        MetadataError::parse(
            path,
            "bgen decompressed probability block length is out of range",
        )
    })?;
    decompressed.clear();
    decompressed.reserve(capacity);
    match compression {
        BgenCompression::Zlib => {
            let mut decoder = ZlibDecoder::new(compressed_payload);
            decoder
                .read_to_end(decompressed)
                .map_err(|source| MetadataError::Io {
                    path: path.to_path_buf(),
                    source,
                })?;
        }
        BgenCompression::Zstd => {
            let decoded = zstd::stream::decode_all(compressed_payload).map_err(|source| {
                MetadataError::Io {
                    path: path.to_path_buf(),
                    source,
                }
            })?;
            decompressed.extend_from_slice(&decoded);
        }
        BgenCompression::None | BgenCompression::Reserved => {
            return Err(MetadataError::parse(
                path,
                "bgen compression value is not a compressed probability block",
            ));
        }
    }

    validate_decompressed_probability_block_len(
        path,
        decompressed.len(),
        expected_decompressed_len,
    )?;
    Ok(())
}

fn validate_decompressed_probability_block_len(
    path: &Path,
    actual_len: usize,
    expected_decompressed_len: u32,
) -> Result<()> {
    let expected_decompressed_len = usize::try_from(expected_decompressed_len).map_err(|_| {
        MetadataError::parse(
            path,
            "bgen decompressed probability block length is out of range",
        )
    })?;
    if actual_len != expected_decompressed_len {
        return Err(MetadataError::parse(
            path,
            "bgen decompressed probability block length does not match length prefix",
        ));
    }
    Ok(())
}

fn skip_exact(reader: &mut impl Read, path: &Path, mut len: u64) -> Result<()> {
    let mut buffer = [0_u8; 8192];
    while len > 0 {
        let chunk_len = buffer
            .len()
            .min(usize::try_from(len).unwrap_or(buffer.len()));
        reader
            .read_exact(&mut buffer[..chunk_len])
            .map_err(|source| MetadataError::Io {
                path: path.to_path_buf(),
                source,
            })?;
        len -= u64::try_from(chunk_len).expect("skip chunk length should fit u64");
    }
    Ok(())
}

fn read_exact_vec(reader: &mut impl Read, path: &Path, len: usize) -> Result<Vec<u8>> {
    let mut bytes = vec![0_u8; len];
    reader
        .read_exact(&mut bytes)
        .map_err(|source| MetadataError::Io {
            path: path.to_path_buf(),
            source,
        })?;
    Ok(bytes)
}

fn read_len_prefixed_string_u16(
    reader: &mut impl Read,
    path: &Path,
    label: &str,
) -> Result<String> {
    let len = usize::from(read_u16_le(reader, path)?);
    read_utf8_string(reader, path, label, len)
}

fn read_len_prefixed_string_u32(
    reader: &mut impl Read,
    path: &Path,
    label: &str,
) -> Result<String> {
    let len = usize::try_from(read_u32_le(reader, path)?)
        .map_err(|_| MetadataError::parse(path, format!("bgen {label} length is out of range")))?;
    read_utf8_string(reader, path, label, len)
}

fn read_utf8_string(
    reader: &mut impl Read,
    path: &Path,
    label: &str,
    len: usize,
) -> Result<String> {
    let mut bytes = vec![0_u8; len];
    reader
        .read_exact(&mut bytes)
        .map_err(|source| MetadataError::Io {
            path: path.to_path_buf(),
            source,
        })?;
    String::from_utf8(bytes)
        .map_err(|error| MetadataError::parse(path, format!("bgen {label} is not UTF-8: {error}")))
}

fn read_u16_le(reader: &mut impl Read, path: &Path) -> Result<u16> {
    let mut bytes = [0_u8; 2];
    reader
        .read_exact(&mut bytes)
        .map_err(|source| MetadataError::Io {
            path: path.to_path_buf(),
            source,
        })?;
    Ok(u16::from_le_bytes(bytes))
}

fn read_u32_le(reader: &mut impl Read, path: &Path) -> Result<u32> {
    let mut bytes = [0_u8; 4];
    reader
        .read_exact(&mut bytes)
        .map_err(|source| MetadataError::Io {
            path: path.to_path_buf(),
            source,
        })?;
    Ok(u32::from_le_bytes(bytes))
}

fn read_u8(reader: &mut impl Read, path: &Path) -> Result<u8> {
    let mut byte = [0_u8; 1];
    reader
        .read_exact(&mut byte)
        .map_err(|source| MetadataError::Io {
            path: path.to_path_buf(),
            source,
        })?;
    Ok(byte[0])
}

#[cfg(test)]
mod tests {
    use super::*;

    fn layout2_payload(
        bit_depth: u8,
        ploidies: &[u8],
        calls: &[Option<(u32, u32)>],
        phased: u8,
    ) -> Vec<u8> {
        assert_eq!(ploidies.len(), calls.len());
        let mut payload = Vec::new();
        payload.extend_from_slice(
            &u32::try_from(calls.len())
                .expect("sample count should fit u32")
                .to_le_bytes(),
        );
        payload.extend_from_slice(&2_u16.to_le_bytes());
        payload.push(2);
        payload.push(2);
        payload.extend_from_slice(ploidies);
        payload.push(phased);
        payload.push(bit_depth);
        append_packed_probabilities(&mut payload, bit_depth, calls);
        payload
    }

    fn append_packed_probabilities(
        output: &mut Vec<u8>,
        bit_depth: u8,
        calls: &[Option<(u32, u32)>],
    ) {
        let mut current_byte = 0_u8;
        let mut bits_in_current_byte = 0_u8;
        for &(p_aa, p_ab) in calls.iter().flatten() {
            append_packed_probability_value(
                output,
                &mut current_byte,
                &mut bits_in_current_byte,
                bit_depth,
                p_aa,
            );
            append_packed_probability_value(
                output,
                &mut current_byte,
                &mut bits_in_current_byte,
                bit_depth,
                p_ab,
            );
        }
        if bits_in_current_byte > 0 {
            output.push(current_byte);
        }
    }

    fn append_packed_probability_value(
        output: &mut Vec<u8>,
        current_byte: &mut u8,
        bits_in_current_byte: &mut u8,
        bit_depth: u8,
        value: u32,
    ) {
        for bit_index in 0..bit_depth {
            let bit = ((value >> bit_index) & 1) as u8;
            *current_byte |= bit << *bits_in_current_byte;
            *bits_in_current_byte += 1;
            if *bits_in_current_byte == 8 {
                output.push(*current_byte);
                *current_byte = 0;
                *bits_in_current_byte = 0;
            }
        }
    }

    fn expected_dosage(bit_depth: u8, p_aa: u32, p_ab: u32) -> f32 {
        let denominator = ((1_u64 << bit_depth) - 1) as f32;
        let p_aa = p_aa as f32 / denominator;
        let p_ab = p_ab as f32 / denominator;
        p_ab + 2.0 * (1.0 - p_aa - p_ab)
    }

    #[test]
    fn decoded_dosage_variant_unpacks_little_endian_probabilities() {
        let payload = layout2_payload(
            3,
            &[2, 0b1000_0010, 2],
            &[Some((5, 2)), None, Some((1, 4))],
            0,
        );
        let decoded = DecodedDosageVariant::decode(Path::new("test.bgen"), &payload, 3, 2)
            .expect("payload should decode");

        let mut values = Vec::new();
        let mut missing = Vec::new();
        decoded
            .decode_source_order(Path::new("test.bgen"), &mut values, &mut missing)
            .expect("dosages should unpack");

        let expected = [expected_dosage(3, 5, 2), 0.0, expected_dosage(3, 1, 4)];
        assert_eq!(missing, vec![false, true, false]);
        for (observed, expected) in values.iter().zip(expected) {
            assert!((observed - expected).abs() <= f32::EPSILON);
        }
    }

    #[test]
    fn decoded_dosage_variant_collapses_phased_haplotype_probabilities() {
        let payload = layout2_payload(
            3,
            &[2, 0b1000_0010, 2],
            &[Some((7, 0)), None, Some((4, 2))],
            1,
        );
        let decoded = DecodedDosageVariant::decode(Path::new("test.bgen"), &payload, 3, 2)
            .expect("payload should decode");

        let mut values = Vec::new();
        let mut missing = Vec::new();
        decoded
            .decode_source_order(Path::new("test.bgen"), &mut values, &mut missing)
            .expect("phased dosages should unpack");

        let denominator = 7.0_f32;
        let expected = [1.0, 0.0, 2.0 - 4.0 / denominator - 2.0 / denominator];
        assert_eq!(missing, vec![false, true, false]);
        for (observed, expected) in values.iter().zip(expected) {
            assert!((observed - expected).abs() <= f32::EPSILON);
        }
    }

    #[test]
    fn decoded_dosage_variant_rejects_impossible_probability_sum() {
        let payload = layout2_payload(1, &[2], &[Some((1, 1))], 0);
        let decoded = DecodedDosageVariant::decode(Path::new("test.bgen"), &payload, 1, 2)
            .expect("payload header should decode");

        let mut values = Vec::new();
        let mut missing = Vec::new();
        let error = decoded
            .decode_source_order(Path::new("test.bgen"), &mut values, &mut missing)
            .expect_err("impossible probabilities should fail");

        assert!(error.to_string().contains("malformed probability"));
    }

    #[test]
    fn decoded_dosage_variant_rejects_trailing_packed_probability_bytes() {
        let mut payload = layout2_payload(8, &[2], &[Some((0, 0))], 0);
        payload.push(0);

        let error = DecodedDosageVariant::decode(Path::new("test.bgen"), &payload, 1, 2)
            .expect_err("trailing packed probability bytes should fail");

        assert!(error.to_string().contains("trailing"));
    }
}
