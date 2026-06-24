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
    pub offsets: Vec<i32>,
    pub values: Vec<u8>,
}

impl StringColumnBuffers {
    pub fn with_capacity(row_capacity: usize, value_capacity: usize) -> Self {
        let mut offsets = Vec::with_capacity(row_capacity.saturating_add(1));
        offsets.push(0);
        Self {
            offsets,
            values: Vec::with_capacity(value_capacity),
        }
    }

    pub fn append_value(&mut self, value: &str) -> Result<(), GenoioError> {
        let next_offset = self.values.len().checked_add(value.len()).ok_or_else(|| {
            GenoioError::unsupported("variant metadata string column exceeds addressable memory")
        })?;
        let next_offset = i32::try_from(next_offset).map_err(|_| {
            GenoioError::unsupported(
                "variant metadata string column exceeds Arrow Utf8 offset capacity",
            )
        })?;
        self.values.extend_from_slice(value.as_bytes());
        self.offsets.push(next_offset);
        Ok(())
    }

    pub fn len(&self) -> usize {
        self.offsets.len().saturating_sub(1)
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

/// Public variant metadata in Arrow-compatible column buffers.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VariantMetadataArrowBuffers {
    pub chroms: StringColumnBuffers,
    pub positions: Vec<i64>,
    pub ids: StringColumnBuffers,
    pub a0s: StringColumnBuffers,
    pub a1s: StringColumnBuffers,
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
        Ok(())
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
    pub samples: Vec<SampleRecord>,
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
