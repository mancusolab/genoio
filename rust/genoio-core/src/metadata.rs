// pattern: Functional Core

use crate::{GenoioError, SourceCapabilities};

/// Sample metadata normalized across VCF, PLINK1, and PLINK2 inputs.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SampleRecord {
    pub fid: Option<String>,
    pub iid: String,
    pub father: Option<String>,
    pub mother: Option<String>,
    pub sex: Option<String>,
    pub phenotype: Option<String>,
    pub source_sample_index: Option<usize>,
    pub haplotype_index: Option<usize>,
}

/// Variant metadata normalized across supported source formats.
#[derive(Debug, Clone, PartialEq)]
pub struct VariantRecord {
    pub chrom: String,
    pub pos: u32,
    pub id: String,
    pub a0: String,
    pub a1: String,
    pub ref_allele: Option<String>,
    pub alt_allele: Option<String>,
    pub source_a0: String,
    pub source_a1: String,
    pub flipped: bool,
    pub qual: Option<f32>,
    pub af: Option<f32>,
    pub maf: Option<f32>,
    pub mac: Option<u32>,
    pub missing_rate: Option<f32>,
    pub n_called: Option<u32>,
}

/// Public sample metadata columns for Python/Arrow adapters.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SampleMetadataColumns {
    pub fids: Vec<Option<String>>,
    pub iids: Vec<String>,
    pub fathers: Vec<Option<String>>,
    pub mothers: Vec<Option<String>>,
    pub sexes: Vec<Option<String>>,
    pub phenotypes: Vec<Option<String>>,
    pub source_sample_indices: Option<Vec<Option<usize>>>,
    pub haplotype_indices: Option<Vec<Option<usize>>>,
}

impl SampleMetadataColumns {
    pub fn from_records(samples: Vec<SampleRecord>, include_haplotype_columns: bool) -> Self {
        let mut fids = Vec::with_capacity(samples.len());
        let mut iids = Vec::with_capacity(samples.len());
        let mut fathers = Vec::with_capacity(samples.len());
        let mut mothers = Vec::with_capacity(samples.len());
        let mut sexes = Vec::with_capacity(samples.len());
        let mut phenotypes = Vec::with_capacity(samples.len());
        let mut source_sample_indices = Vec::with_capacity(samples.len());
        let mut haplotype_indices = Vec::with_capacity(samples.len());

        for sample in samples {
            fids.push(sample.fid);
            iids.push(sample.iid);
            fathers.push(sample.father);
            mothers.push(sample.mother);
            sexes.push(sample.sex);
            phenotypes.push(sample.phenotype);
            if include_haplotype_columns {
                source_sample_indices.push(sample.source_sample_index);
                haplotype_indices.push(sample.haplotype_index);
            }
        }

        Self {
            fids,
            iids,
            fathers,
            mothers,
            sexes,
            phenotypes,
            source_sample_indices: include_haplotype_columns.then_some(source_sample_indices),
            haplotype_indices: include_haplotype_columns.then_some(haplotype_indices),
        }
    }
}

/// Arrow-compatible nullable UTF-8 column buffers for Python/Arrow adapters.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NullableStringColumnBuffers {
    /// Arrow Utf8 byte offsets with a leading zero and one offset per row end.
    pub offsets: Vec<i32>,
    /// Contiguous UTF-8 bytes for all non-null values.
    pub values: Vec<u8>,
    /// Per-row validity bits; false rows reuse the previous byte offset.
    pub validity: Vec<bool>,
}

impl NullableStringColumnBuffers {
    /// Allocate a nullable string column with Arrow-compatible initial offset state.
    pub fn with_capacity(row_capacity: usize, value_capacity: usize) -> Self {
        let mut offsets = Vec::with_capacity(row_capacity.saturating_add(1));
        offsets.push(0);
        Self {
            offsets,
            values: Vec::with_capacity(value_capacity),
            validity: Vec::with_capacity(row_capacity),
        }
    }

    /// Append one optional UTF-8 value while preserving Arrow offset invariants.
    pub fn append_option(&mut self, value: Option<&str>) -> Result<(), GenoioError> {
        match value {
            Some(value) => {
                append_utf8_value(&mut self.offsets, &mut self.values, value)?;
                self.validity.push(true);
            }
            None => {
                self.offsets.push(*self.offsets.last().ok_or_else(|| {
                    GenoioError::internal_contract("nullable string column missing initial offset")
                })?);
                self.validity.push(false);
            }
        }
        Ok(())
    }

    /// Return the number of logical rows in the column.
    pub fn len(&self) -> usize {
        self.validity.len()
    }

    /// Return true when the column contains no logical rows.
    pub fn is_empty(&self) -> bool {
        self.validity.is_empty()
    }
}

/// Sample metadata with public fields staged in Arrow-compatible column buffers.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SampleMetadataArrowBuffers {
    /// Family IDs; nullable because VCF/BGEN sources may not provide them.
    pub fids: NullableStringColumnBuffers,
    /// Individual/sample IDs; required for every public sample row.
    pub iids: StringColumnBuffers,
    /// Paternal IDs from pedigree-style metadata, when present.
    pub fathers: NullableStringColumnBuffers,
    /// Maternal IDs from pedigree-style metadata, when present.
    pub mothers: NullableStringColumnBuffers,
    /// Source sex labels, when present.
    pub sexes: NullableStringColumnBuffers,
    /// Source phenotype labels, when present.
    pub phenotypes: NullableStringColumnBuffers,
    /// Source sample row for haplotype-expanded outputs; absent for genotype outputs.
    pub source_sample_indices: Option<Vec<Option<usize>>>,
    /// Haplotype row index for haplotype-expanded outputs; absent for genotype outputs.
    pub haplotype_indices: Option<Vec<Option<usize>>>,
}

/// Borrowed view of a string inside columnar sample metadata.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SampleMetadataStr<'a>(&'a str);

impl<'a> SampleMetadataStr<'a> {
    /// Return the underlying UTF-8 value.
    pub fn as_str(self) -> &'a str {
        self.0
    }
}

/// Borrowed row view produced from columnar sample metadata on demand.
///
/// This exists for tests and diagnostics that need row-shaped assertions. The
/// main reader and Python paths should keep using `SampleMetadataArrowBuffers`
/// directly so they do not re-materialize sample records.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SampleMetadataRow<'a> {
    pub fid: Option<SampleMetadataStr<'a>>,
    pub iid: SampleMetadataStr<'a>,
    pub father: Option<SampleMetadataStr<'a>>,
    pub mother: Option<SampleMetadataStr<'a>>,
    pub sex: Option<SampleMetadataStr<'a>>,
    pub phenotype: Option<SampleMetadataStr<'a>>,
    pub source_sample_index: Option<usize>,
    pub haplotype_index: Option<usize>,
}

#[derive(Debug, Clone)]
pub struct SampleMetadataRows<'a> {
    buffers: &'a SampleMetadataArrowBuffers,
    index: usize,
}

impl SampleMetadataArrowBuffers {
    /// Allocate sample metadata columns for `row_capacity` public sample rows.
    ///
    /// `include_haplotype_columns` controls whether the public Arrow schema
    /// includes source-sample and haplotype mapping columns. Genotype outputs
    /// leave those columns absent to preserve the existing public schema.
    pub fn with_capacity(row_capacity: usize, include_haplotype_columns: bool) -> Self {
        Self {
            fids: NullableStringColumnBuffers::with_capacity(
                row_capacity,
                row_capacity.saturating_mul(8),
            ),
            iids: StringColumnBuffers::with_capacity(row_capacity, row_capacity.saturating_mul(12)),
            fathers: NullableStringColumnBuffers::with_capacity(
                row_capacity,
                row_capacity.saturating_mul(8),
            ),
            mothers: NullableStringColumnBuffers::with_capacity(
                row_capacity,
                row_capacity.saturating_mul(8),
            ),
            sexes: NullableStringColumnBuffers::with_capacity(
                row_capacity,
                row_capacity.saturating_mul(8),
            ),
            phenotypes: NullableStringColumnBuffers::with_capacity(
                row_capacity,
                row_capacity.saturating_mul(8),
            ),
            source_sample_indices: include_haplotype_columns
                .then(|| Vec::with_capacity(row_capacity)),
            haplotype_indices: include_haplotype_columns.then(|| Vec::with_capacity(row_capacity)),
        }
    }

    /// Append one normalized sample record into the columnar buffers.
    pub fn push_record(&mut self, sample: &SampleRecord) -> Result<(), GenoioError> {
        self.fids.append_option(sample.fid.as_deref())?;
        self.iids.append_value(&sample.iid)?;
        self.fathers.append_option(sample.father.as_deref())?;
        self.mothers.append_option(sample.mother.as_deref())?;
        self.sexes.append_option(sample.sex.as_deref())?;
        self.phenotypes.append_option(sample.phenotype.as_deref())?;
        if let Some(source_sample_indices) = self.source_sample_indices.as_mut() {
            source_sample_indices.push(sample.source_sample_index);
        }
        if let Some(haplotype_indices) = self.haplotype_indices.as_mut() {
            haplotype_indices.push(sample.haplotype_index);
        }
        Ok(())
    }

    /// Build Arrow-compatible public sample buffers from normalized records.
    pub fn from_records(
        samples: &[SampleRecord],
        include_haplotype_columns: bool,
    ) -> Result<Self, GenoioError> {
        let mut output = Self::with_capacity(samples.len(), include_haplotype_columns);
        for sample in samples {
            output.push_record(sample)?;
        }
        Ok(output)
    }

    /// Build sample buffers only when the caller requested sample metadata.
    ///
    /// Returning `None` distinguishes omitted metadata from a requested but
    /// empty sample table, which is important at the PyO3 boundary.
    pub fn optional_from_records(
        samples: &[SampleRecord],
        return_samples: bool,
        include_haplotype_columns: bool,
    ) -> Result<Option<Self>, GenoioError> {
        if !return_samples {
            return Ok(None);
        }
        Self::from_records(samples, include_haplotype_columns).map(Some)
    }

    /// Return the number of public sample rows represented by these buffers.
    pub fn len(&self) -> usize {
        self.iids.len()
    }

    /// Return true when no public sample rows are present.
    pub fn is_empty(&self) -> bool {
        self.iids.is_empty()
    }

    /// Iterate over borrowed row views without changing the columnar storage.
    pub fn iter(&self) -> SampleMetadataRows<'_> {
        SampleMetadataRows {
            buffers: self,
            index: 0,
        }
    }
}

/// Public variant metadata columns for Python/Arrow adapters.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VariantMetadataColumns {
    pub chroms: Vec<String>,
    pub positions: Vec<u32>,
    pub ids: Vec<String>,
    pub a0s: Vec<String>,
    pub a1s: Vec<String>,
}

impl VariantMetadataColumns {
    pub fn from_records(variants: Vec<VariantRecord>) -> Self {
        let mut chroms = Vec::with_capacity(variants.len());
        let mut positions = Vec::with_capacity(variants.len());
        let mut ids = Vec::with_capacity(variants.len());
        let mut a0s = Vec::with_capacity(variants.len());
        let mut a1s = Vec::with_capacity(variants.len());

        for variant in variants {
            chroms.push(variant.chrom);
            positions.push(variant.pos);
            ids.push(variant.id);
            a0s.push(variant.a0);
            a1s.push(variant.a1);
        }

        Self {
            chroms,
            positions,
            ids,
            a0s,
            a1s,
        }
    }
}

/// Arrow-compatible UTF-8 column buffers for Python/Arrow adapters.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StringColumnBuffers {
    /// Arrow Utf8 byte offsets with a leading zero and one offset per row end.
    pub offsets: Vec<i32>,
    /// Contiguous UTF-8 bytes for all values.
    pub values: Vec<u8>,
}

impl StringColumnBuffers {
    /// Allocate a non-null string column with Arrow-compatible initial offset state.
    pub fn with_capacity(row_capacity: usize, value_capacity: usize) -> Self {
        let mut offsets = Vec::with_capacity(row_capacity.saturating_add(1));
        offsets.push(0);
        Self {
            offsets,
            values: Vec::with_capacity(value_capacity),
        }
    }

    /// Append one non-null UTF-8 value while preserving Arrow offset invariants.
    pub fn append_value(&mut self, value: &str) -> Result<(), GenoioError> {
        append_utf8_value(&mut self.offsets, &mut self.values, value)
    }

    /// Return the number of logical rows in the column.
    pub fn len(&self) -> usize {
        self.offsets.len().saturating_sub(1)
    }

    /// Return true when the column contains no logical rows.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

impl<'a> Iterator for SampleMetadataRows<'a> {
    type Item = SampleMetadataRow<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.index >= self.buffers.len() {
            return None;
        }
        let index = self.index;
        self.index += 1;
        Some(SampleMetadataRow {
            fid: nullable_string_at(&self.buffers.fids, index),
            iid: SampleMetadataStr(string_at(&self.buffers.iids, index)),
            father: nullable_string_at(&self.buffers.fathers, index),
            mother: nullable_string_at(&self.buffers.mothers, index),
            sex: nullable_string_at(&self.buffers.sexes, index),
            phenotype: nullable_string_at(&self.buffers.phenotypes, index),
            source_sample_index: optional_usize_at(&self.buffers.source_sample_indices, index),
            haplotype_index: optional_usize_at(&self.buffers.haplotype_indices, index),
        })
    }
}

fn string_at(column: &StringColumnBuffers, index: usize) -> &str {
    let start = column.offsets[index] as usize;
    let end = column.offsets[index + 1] as usize;
    // SAFETY: StringColumnBuffers is populated only from `&str` values, and
    // offsets are appended immediately after those valid UTF-8 bytes.
    unsafe { std::str::from_utf8_unchecked(&column.values[start..end]) }
}

fn nullable_string_at(
    column: &NullableStringColumnBuffers,
    index: usize,
) -> Option<SampleMetadataStr<'_>> {
    column.validity[index].then(|| {
        let start = column.offsets[index] as usize;
        let end = column.offsets[index + 1] as usize;
        // SAFETY: NullableStringColumnBuffers is populated only from `&str`
        // values, and offsets are appended immediately after valid UTF-8 bytes.
        SampleMetadataStr(unsafe { std::str::from_utf8_unchecked(&column.values[start..end]) })
    })
}

fn optional_usize_at(values: &Option<Vec<Option<usize>>>, index: usize) -> Option<usize> {
    values.as_ref().and_then(|values| values[index])
}

fn append_utf8_value(
    offsets: &mut Vec<i32>,
    values: &mut Vec<u8>,
    value: &str,
) -> Result<(), GenoioError> {
    // Store values in Arrow Utf8 layout now so PyO3 can hand ownership to
    // Arrow arrays without rebuilding strings row by row.
    let next_offset = values.len().checked_add(value.len()).ok_or_else(|| {
        GenoioError::unsupported("metadata string column exceeds addressable memory")
    })?;
    let next_offset = i32::try_from(next_offset).map_err(|_| {
        GenoioError::unsupported("metadata string column exceeds Arrow Utf8 offset capacity")
    })?;
    values.extend_from_slice(value.as_bytes());
    offsets.push(next_offset);
    Ok(())
}

/// Variant metadata with public fields staged in Arrow-compatible column buffers.
///
/// Non-public fields are retained as columns so format readers can attach
/// genotype-derived statistics without leaving the columnar backend.
#[derive(Debug, Clone, PartialEq)]
pub struct VariantMetadataArrowBuffers {
    pub chroms: StringColumnBuffers,
    pub positions: Vec<i64>,
    pub ids: StringColumnBuffers,
    pub a0s: StringColumnBuffers,
    pub a1s: StringColumnBuffers,
    pub ref_alleles: Vec<Option<String>>,
    pub alt_alleles: Vec<Option<String>>,
    pub source_a0s: StringColumnBuffers,
    pub source_a1s: StringColumnBuffers,
    pub flipped: Vec<bool>,
    pub quals: Vec<Option<f32>>,
    pub afs: Vec<Option<f32>>,
    pub mafs: Vec<Option<f32>>,
    pub macs: Vec<Option<u32>>,
    pub missing_rates: Vec<Option<f32>>,
    pub n_called: Vec<Option<u32>>,
}

impl VariantMetadataArrowBuffers {
    pub fn with_capacity(row_capacity: usize) -> Self {
        Self {
            chroms: StringColumnBuffers::with_capacity(
                row_capacity,
                row_capacity.saturating_mul(8),
            ),
            positions: Vec::with_capacity(row_capacity),
            ids: StringColumnBuffers::with_capacity(row_capacity, row_capacity.saturating_mul(12)),
            a0s: StringColumnBuffers::with_capacity(row_capacity, row_capacity.saturating_mul(2)),
            a1s: StringColumnBuffers::with_capacity(row_capacity, row_capacity.saturating_mul(2)),
            ref_alleles: Vec::with_capacity(row_capacity),
            alt_alleles: Vec::with_capacity(row_capacity),
            source_a0s: StringColumnBuffers::with_capacity(
                row_capacity,
                row_capacity.saturating_mul(2),
            ),
            source_a1s: StringColumnBuffers::with_capacity(
                row_capacity,
                row_capacity.saturating_mul(2),
            ),
            flipped: Vec::with_capacity(row_capacity),
            quals: Vec::with_capacity(row_capacity),
            afs: Vec::with_capacity(row_capacity),
            mafs: Vec::with_capacity(row_capacity),
            macs: Vec::with_capacity(row_capacity),
            missing_rates: Vec::with_capacity(row_capacity),
            n_called: Vec::with_capacity(row_capacity),
        }
    }

    pub fn push(
        &mut self,
        chrom: &str,
        pos: i64,
        id: &str,
        a0: &str,
        a1: &str,
    ) -> Result<(), GenoioError> {
        self.chroms.append_value(chrom)?;
        self.positions.push(pos);
        self.ids.append_value(id)?;
        self.a0s.append_value(a0)?;
        self.a1s.append_value(a1)?;
        self.ref_alleles.push(None);
        self.alt_alleles.push(None);
        self.source_a0s.append_value(a0)?;
        self.source_a1s.append_value(a1)?;
        self.flipped.push(false);
        self.quals.push(None);
        self.afs.push(None);
        self.mafs.push(None);
        self.macs.push(None);
        self.missing_rates.push(None);
        self.n_called.push(None);
        Ok(())
    }

    pub fn push_record(&mut self, variant: &VariantRecord) -> Result<(), GenoioError> {
        self.chroms.append_value(&variant.chrom)?;
        self.positions.push(i64::from(variant.pos));
        self.ids.append_value(&variant.id)?;
        self.a0s.append_value(&variant.a0)?;
        self.a1s.append_value(&variant.a1)?;
        self.ref_alleles.push(variant.ref_allele.clone());
        self.alt_alleles.push(variant.alt_allele.clone());
        self.source_a0s.append_value(&variant.source_a0)?;
        self.source_a1s.append_value(&variant.source_a1)?;
        self.flipped.push(variant.flipped);
        self.quals.push(variant.qual);
        self.afs.push(variant.af);
        self.mafs.push(variant.maf);
        self.macs.push(variant.mac);
        self.missing_rates.push(variant.missing_rate);
        self.n_called.push(variant.n_called);
        Ok(())
    }

    /// Build Arrow-compatible public variant buffers from normalized records.
    pub fn from_records(variants: &[VariantRecord]) -> Result<Self, GenoioError> {
        let mut output = Self::with_capacity(variants.len());
        for variant in variants {
            output.push_record(variant)?;
        }
        Ok(output)
    }

    pub fn len(&self) -> usize {
        self.positions.len()
    }

    pub fn is_empty(&self) -> bool {
        self.positions.is_empty()
    }
}

/// Complete source metadata with public variants already staged as column buffers.
#[derive(Debug, Clone, PartialEq)]
pub struct MetadataArrowOutput {
    pub samples: SampleMetadataArrowBuffers,
    pub variants: VariantMetadataArrowBuffers,
    pub capabilities: SourceCapabilities,
}

/// Complete source metadata returned before or alongside matrix reads.
#[derive(Debug, Clone, PartialEq)]
pub struct MetadataOutput {
    pub samples: Vec<SampleRecord>,
    pub variants: Vec<VariantRecord>,
    pub capabilities: SourceCapabilities,
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample(iid: &str, source_sample_index: Option<usize>) -> SampleRecord {
        SampleRecord {
            fid: Some(format!("F_{iid}")),
            iid: iid.to_string(),
            father: None,
            mother: Some("0".to_string()),
            sex: Some("unknown".to_string()),
            phenotype: None,
            source_sample_index,
            haplotype_index: None,
        }
    }

    fn variant(id: &str, pos: u32) -> VariantRecord {
        VariantRecord {
            chrom: "1".to_string(),
            pos,
            id: id.to_string(),
            a0: "A".to_string(),
            a1: "G".to_string(),
            ref_allele: Some("A".to_string()),
            alt_allele: Some("G,T".to_string()),
            source_a0: "A".to_string(),
            source_a1: "G".to_string(),
            flipped: false,
            qual: Some(99.0),
            af: Some(0.2),
            maf: Some(0.2),
            mac: Some(4),
            missing_rate: Some(0.0),
            n_called: Some(20),
        }
    }

    #[test]
    fn sample_metadata_columns_preserve_public_order_and_optional_mapping_columns() {
        let columns = SampleMetadataColumns::from_records(
            vec![sample("s1", Some(2)), sample("s2", Some(7))],
            true,
        );

        assert_eq!(
            columns.fids,
            vec![Some("F_s1".to_string()), Some("F_s2".to_string())]
        );
        assert_eq!(columns.iids, vec!["s1".to_string(), "s2".to_string()]);
        assert_eq!(
            columns.mothers,
            vec![Some("0".to_string()), Some("0".to_string())]
        );
        assert_eq!(columns.source_sample_indices, Some(vec![Some(2), Some(7)]));
        assert_eq!(columns.haplotype_indices, Some(vec![None, None]));

        let columns_without_mapping =
            SampleMetadataColumns::from_records(vec![sample("s3", Some(9))], false);
        assert_eq!(columns_without_mapping.source_sample_indices, None);
        assert_eq!(columns_without_mapping.haplotype_indices, None);
    }

    #[test]
    fn sample_metadata_arrow_buffers_preserve_nullable_strings_and_mapping_columns() {
        let buffers = SampleMetadataArrowBuffers::from_records(
            &[sample("s1", Some(2)), sample("s2", Some(7))],
            true,
        )
        .expect("sample buffers should be built");

        assert_eq!(buffers.len(), 2);
        assert_eq!(buffers.fids.validity, vec![true, true]);
        assert_eq!(buffers.iids.len(), 2);
        assert_eq!(buffers.mothers.validity, vec![true, true]);
        assert_eq!(buffers.phenotypes.validity, vec![false, false]);
        assert_eq!(buffers.source_sample_indices, Some(vec![Some(2), Some(7)]));
        assert_eq!(buffers.haplotype_indices, Some(vec![None, None]));

        let buffers_without_mapping =
            SampleMetadataArrowBuffers::from_records(&[sample("s3", Some(9))], false)
                .expect("sample buffers should be built");
        assert_eq!(buffers_without_mapping.source_sample_indices, None);
        assert_eq!(buffers_without_mapping.haplotype_indices, None);
    }

    #[test]
    fn variant_metadata_columns_preserve_public_order_and_drop_internal_fields() {
        let columns =
            VariantMetadataColumns::from_records(vec![variant("rs1", 10), variant("rs2", 20)]);

        assert_eq!(columns.chroms, vec!["1".to_string(), "1".to_string()]);
        assert_eq!(columns.positions, vec![10, 20]);
        assert_eq!(columns.ids, vec!["rs1".to_string(), "rs2".to_string()]);
        assert_eq!(columns.a0s, vec!["A".to_string(), "A".to_string()]);
        assert_eq!(columns.a1s, vec!["G".to_string(), "G".to_string()]);
    }

    #[test]
    fn string_column_buffers_append_arrow_utf8_offsets_and_values() {
        let mut buffers = StringColumnBuffers::with_capacity(2, 8);

        buffers.append_value("chr1").unwrap();
        buffers.append_value("A").unwrap();

        assert_eq!(buffers.offsets, vec![0, 4, 5]);
        assert_eq!(buffers.values, b"chr1A");
        assert_eq!(buffers.len(), 2);
    }

    #[test]
    fn variant_metadata_arrow_buffers_preserve_public_columns() {
        let mut buffers = VariantMetadataArrowBuffers::with_capacity(2);

        buffers.push("1", 10, "rs1", "A", "G").unwrap();
        buffers.push("2", 20, ".", "C", "T").unwrap();

        assert_eq!(buffers.chroms.offsets, vec![0, 1, 2]);
        assert_eq!(buffers.chroms.values, b"12");
        assert_eq!(buffers.positions, vec![10, 20]);
        assert_eq!(buffers.ids.offsets, vec![0, 3, 4]);
        assert_eq!(buffers.ids.values, b"rs1.");
        assert_eq!(buffers.a0s.values, b"AC");
        assert_eq!(buffers.a1s.values, b"GT");
    }
}
