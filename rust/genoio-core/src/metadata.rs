// pattern: Functional Core

use crate::{
    filter::{VariantMetadataView, VariantStats},
    GenoioError, SourceCapabilities,
};

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
pub struct SampleMetadataBuffers {
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
/// main reader and Python paths should keep using `SampleMetadataBuffers`
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
    buffers: &'a SampleMetadataBuffers,
    index: usize,
}

impl SampleMetadataBuffers {
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

    /// Replace one logical row while preserving Arrow Utf8 offsets.
    ///
    /// Retained-row allele flips can swap strings after a row has already been
    /// appended. This keeps the column in canonical Arrow Utf8 layout instead
    /// of rebuilding the entire metadata buffer.
    pub fn replace_value(&mut self, index: usize, value: &str) -> Result<(), GenoioError> {
        if index >= self.len() {
            return Err(metadata_row_index_error(index, self.len()));
        }

        let start = self.offsets[index] as usize;
        let end = self.offsets[index + 1] as usize;
        let old_len = end.saturating_sub(start);
        let new_len = value.len();
        if old_len == new_len {
            self.values[start..end].copy_from_slice(value.as_bytes());
            return Ok(());
        }

        let new_total_len = self
            .values
            .len()
            .checked_sub(old_len)
            .and_then(|len| len.checked_add(new_len))
            .ok_or_else(|| {
                GenoioError::unsupported("metadata string column exceeds addressable memory")
            })?;
        i32::try_from(new_total_len).map_err(|_| {
            GenoioError::unsupported("metadata string column exceeds Arrow Utf8 offset capacity")
        })?;

        // Allele strings may have different byte lengths, so every following
        // offset must move by the replacement delta.
        self.values.splice(start..end, value.bytes());
        let delta = i64::try_from(new_len)
            .and_then(|new_len| i64::try_from(old_len).map(|old_len| new_len - old_len))
            .map_err(|_| {
                GenoioError::unsupported("metadata string column exceeds addressable memory")
            })?;
        for offset in &mut self.offsets[index + 1..] {
            let next_offset = i64::from(*offset) + delta;
            *offset = i32::try_from(next_offset).map_err(|_| {
                GenoioError::unsupported(
                    "metadata string column exceeds Arrow Utf8 offset capacity",
                )
            })?;
        }
        Ok(())
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
pub struct VariantMetadataBuffers {
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

impl VariantMetadataBuffers {
    /// Allocate variant metadata columns for `row_capacity` public variant rows.
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

    /// Append one owned compatibility record by borrowing its metadata fields.
    pub fn push_record(&mut self, variant: &VariantRecord) -> Result<(), GenoioError> {
        self.push_view(variant)
    }

    /// Append one borrowed metadata view into Arrow-compatible columns.
    ///
    /// This is the columnar hot-path boundary used by parser-specific phases:
    /// callers can project borrowed parser fields into Arrow buffers without
    /// first materializing an owned [`VariantRecord`].
    pub fn push_view<V: VariantMetadataView + ?Sized>(
        &mut self,
        variant: &V,
    ) -> Result<(), GenoioError> {
        self.chroms.append_value(variant.chrom())?;
        self.positions.push(i64::from(variant.pos()));
        self.ids.append_value(variant.id())?;
        self.a0s.append_value(variant.a0())?;
        self.a1s.append_value(variant.a1())?;
        self.ref_alleles
            .push(variant.ref_allele().map(str::to_owned));
        self.alt_alleles
            .push(variant.alt_allele().map(str::to_owned));
        self.source_a0s.append_value(variant.source_a0())?;
        self.source_a1s.append_value(variant.source_a1())?;
        self.flipped.push(variant.flipped());
        self.quals.push(variant.qual());
        self.afs.push(variant.af());
        self.mafs.push(variant.maf());
        self.macs.push(variant.mac());
        self.missing_rates.push(variant.missing_rate());
        self.n_called.push(variant.n_called());
        Ok(())
    }

    /// Append one borrowed metadata view and attach genotype statistics.
    ///
    /// Use this when a genotype-stat filter has already retained the row and
    /// the public variant metadata should expose the computed statistics.
    pub fn push_view_with_stats<V: VariantMetadataView + ?Sized>(
        &mut self,
        variant: &V,
        stats: VariantStats,
    ) -> Result<(), GenoioError> {
        self.push_view(variant)?;
        self.attach_stats(self.len() - 1, stats)
    }

    /// Attach computed genotype statistics to a retained metadata row.
    ///
    /// The row index is in retained-output coordinates, not source-row
    /// coordinates. This keeps delayed stat attachment aligned with the matrix
    /// column that survived filtering and windowing.
    pub fn attach_stats(
        &mut self,
        row_index: usize,
        stats: VariantStats,
    ) -> Result<(), GenoioError> {
        self.validate_row_index(row_index)?;
        self.afs[row_index] = stats.af.map(|value| value as f32);
        self.mafs[row_index] = stats.maf.map(|value| value as f32);
        // Public variant metadata keeps MAC integer-valued. Dosage filters can
        // still evaluate fractional MAC internally before this projection.
        self.macs[row_index] = stats.mac.and_then(exact_u32_from_f64);
        self.missing_rates[row_index] = Some(stats.missing_rate as f32);
        self.n_called[row_index] = Some(stats.n_called);
        Ok(())
    }

    /// Flip public alleles for a retained row when values are encoded as minor allele.
    ///
    /// This mutates only the public allele columns and flip marker. Source
    /// allele columns remain unchanged so downstream adapters can report both
    /// public and original source orientation.
    pub fn flip_to_minor_allele(&mut self, row_index: usize) -> Result<(), GenoioError> {
        self.validate_row_index(row_index)?;
        let a0 = string_at(&self.a0s, row_index).to_owned();
        let a1 = string_at(&self.a1s, row_index).to_owned();
        self.a0s.replace_value(row_index, &a1)?;
        self.a1s.replace_value(row_index, &a0)?;
        self.flipped[row_index] = true;
        if let Some(af) = self.afs[row_index] {
            self.afs[row_index] = Some(1.0 - af);
        }
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

    fn validate_row_index(&self, row_index: usize) -> Result<(), GenoioError> {
        if row_index < self.len() {
            Ok(())
        } else {
            Err(metadata_row_index_error(row_index, self.len()))
        }
    }
}

fn metadata_row_index_error(index: usize, len: usize) -> GenoioError {
    GenoioError::internal_contract(format!(
        "metadata row index {index} is outside row count {len}"
    ))
}

fn exact_u32_from_f64(value: f64) -> Option<u32> {
    if value.is_finite() && value.fract() == 0.0 && value >= 0.0 && value <= f64::from(u32::MAX) {
        Some(value as u32)
    } else {
        None
    }
}

/// Complete source metadata with public variants already staged as column buffers.
#[derive(Debug, Clone, PartialEq)]
pub struct MetadataOutput {
    pub samples: SampleMetadataBuffers,
    pub variants: VariantMetadataBuffers,
    pub capabilities: SourceCapabilities,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::filter::{VariantMetadataView, VariantStats};

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

    struct BorrowedVariantView<'a> {
        chrom: &'a str,
        pos: u32,
        id: &'a str,
        a0: &'a str,
        a1: &'a str,
        ref_allele: Option<&'a str>,
        alt_allele: Option<&'a str>,
        source_a0: &'a str,
        source_a1: &'a str,
        qual: Option<f32>,
    }

    impl VariantMetadataView for BorrowedVariantView<'_> {
        fn chrom(&self) -> &str {
            self.chrom
        }

        fn pos(&self) -> u32 {
            self.pos
        }

        fn id(&self) -> &str {
            self.id
        }

        fn a0(&self) -> &str {
            self.a0
        }

        fn a1(&self) -> &str {
            self.a1
        }

        fn ref_allele(&self) -> Option<&str> {
            self.ref_allele
        }

        fn alt_allele(&self) -> Option<&str> {
            self.alt_allele
        }

        fn source_a0(&self) -> &str {
            self.source_a0
        }

        fn source_a1(&self) -> &str {
            self.source_a1
        }

        fn qual(&self) -> Option<f32> {
            self.qual
        }
    }

    fn borrowed_variant() -> BorrowedVariantView<'static> {
        BorrowedVariantView {
            chrom: "chr1",
            pos: 42,
            id: "rs42",
            a0: "A",
            a1: "G",
            ref_allele: Some("A"),
            alt_allele: Some("G"),
            source_a0: "A",
            source_a1: "G",
            qual: Some(30.0),
        }
    }

    #[test]
    fn sample_metadata_buffers_preserve_nullable_strings_and_mapping_columns() {
        let buffers = SampleMetadataBuffers::from_records(
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
            SampleMetadataBuffers::from_records(&[sample("s3", Some(9))], false)
                .expect("sample buffers should be built");
        assert_eq!(buffers_without_mapping.source_sample_indices, None);
        assert_eq!(buffers_without_mapping.haplotype_indices, None);
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
    fn variant_metadata_buffers_preserve_public_columns() {
        let mut buffers = VariantMetadataBuffers::with_capacity(2);

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

    #[test]
    fn variant_metadata_buffers_append_borrowed_views() {
        let mut buffers = VariantMetadataBuffers::with_capacity(1);

        buffers.push_view(&borrowed_variant()).unwrap();

        assert_eq!(buffers.positions, vec![42]);
        assert_eq!(buffers.ref_alleles, vec![Some("A".to_string())]);
        assert_eq!(buffers.alt_alleles, vec![Some("G".to_string())]);
        assert_eq!(buffers.quals, vec![Some(30.0)]);
        assert_eq!(buffers.afs, vec![None]);
    }

    #[test]
    fn variant_metadata_buffers_append_views_with_stats() {
        let mut buffers = VariantMetadataBuffers::with_capacity(1);
        let stats = VariantStats {
            af: Some(0.75),
            maf: Some(0.25),
            mac: Some(3.0),
            missing_rate: 0.125,
            n_called: 12,
            polymorphic: true,
        };

        buffers
            .push_view_with_stats(&borrowed_variant(), stats)
            .unwrap();

        assert_eq!(buffers.afs, vec![Some(0.75)]);
        assert_eq!(buffers.mafs, vec![Some(0.25)]);
        assert_eq!(buffers.macs, vec![Some(3)]);
        assert_eq!(buffers.missing_rates, vec![Some(0.125)]);
        assert_eq!(buffers.n_called, vec![Some(12)]);
    }

    #[test]
    fn variant_metadata_buffers_mutate_retained_rows() {
        let mut buffers = VariantMetadataBuffers::with_capacity(1);
        buffers.push_view(&borrowed_variant()).unwrap();

        buffers.flip_to_minor_allele(0).unwrap();
        buffers
            .attach_stats(
                0,
                VariantStats {
                    af: Some(0.25),
                    maf: Some(0.25),
                    mac: Some(1.5),
                    missing_rate: 0.0,
                    n_called: 4,
                    polymorphic: true,
                },
            )
            .unwrap();

        assert_eq!(string_at(&buffers.a0s, 0), "G");
        assert_eq!(string_at(&buffers.a1s, 0), "A");
        assert_eq!(buffers.flipped, vec![true]);
        assert_eq!(buffers.afs, vec![Some(0.25)]);
        assert_eq!(buffers.macs, vec![None]);
    }
}
