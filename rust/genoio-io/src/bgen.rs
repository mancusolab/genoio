// pattern: Imperative Shell
//! BGEN reader orchestration and matrix assembly.

use std::fs::File;
use std::io::{BufReader, Seek, SeekFrom};
use std::path::Path;

use crate::matrix::{
    finish_dense_matrix, finish_variant_major_dense_matrix, shrink_sample_major_width,
    write_sample_major_variant_slot, DenseMatrixParts, VariantMajorDenseParts,
};
use crate::retention::{MetadataRetentionAction, RetainedVariantState, RetentionAction};
use crate::Result;
use genoio_core::{
    attach_variant_stats, compute_dosage_variant_stats, is_dosage_polymorphic,
    select_samples_source_order, DenseGenotypeMatrix, DenseSampleSelection, GenoioError,
    GenotypeFilterConjunction, GenotypeFilterPlan, MetadataOutput, PartialFilterDecision,
    SampleRecord, SourceCapabilities, VariantFilter, VariantRecord, VariantStats, VariantWindow,
};

mod decode;
mod header;
mod index;
mod io;

const BGEN_READER_BUFFER_SIZE: usize = 1 << 20;

use decode::{
    decode_buffered_dosage_values, decode_buffered_haplotype_values,
    read_layout2_probability_payload_into, skip_layout2_probability_payload,
    try_decode_buffered_dosage_values_into_sample_major_slot, DosageDecodeBuffers,
    HaplotypeDecodeBuffers, SampleMajorSlotMut,
};
use header::{
    read_bgen_samples, read_layout2_variant_identifying_data, read_layout2_variant_metadata,
    skip_layout2_variant_identifying_data, BgenHeader,
};
use index::{indexed_region_records, validate_index_record_consumed, BgenIndexRecord};

#[derive(Debug, Clone, Copy, Default)]
struct DosageFilterCounts {
    allele_count: f64,
    called_count: u64,
    missing_count: u64,
}

impl DosageFilterCounts {
    fn evaluate_plan(self, plan: GenotypeFilterPlan) -> Result<Option<bool>> {
        match plan {
            GenotypeFilterPlan::Generic => Ok(None),
            GenotypeFilterPlan::Polymorphic => Ok(Some(self.is_polymorphic()?)),
            GenotypeFilterPlan::MacRange { min, max } => Ok(Some(self.mac_in_range(min, max)?)),
            GenotypeFilterPlan::MafRange { min, max } => Ok(Some(self.maf_in_range(min, max)?)),
            GenotypeFilterPlan::MissingRateMax { max } => {
                Ok(Some(self.missing_rate()? <= f64::from(max)))
            }
            GenotypeFilterPlan::Conjunction(plan) => Ok(Some(self.evaluate_conjunction(plan)?)),
        }
    }

    fn evaluate_conjunction(self, plan: GenotypeFilterConjunction) -> Result<bool> {
        if plan.polymorphic && !self.is_polymorphic()? {
            return Ok(false);
        }
        if (plan.mac_min.is_some() || plan.mac_max.is_some())
            && !self.mac_in_range(plan.mac_min, plan.mac_max)?
        {
            return Ok(false);
        }
        if (plan.maf_min.is_some() || plan.maf_max.is_some())
            && !self.maf_in_range(plan.maf_min, plan.maf_max)?
        {
            return Ok(false);
        }
        if let Some(max) = plan.missing_rate_max {
            if self.missing_rate()? > f64::from(max) {
                return Ok(false);
            }
        }
        Ok(true)
    }

    fn called_alleles(self) -> Result<Option<f64>> {
        if self.called_count == 0 {
            return Ok(None);
        }
        let called_count = u32::try_from(self.called_count).map_err(|_| {
            GenoioError::invalid_source(
                "<filter>",
                "called genotype count exceeds supported metadata range",
            )
        })?;
        Ok(Some(2.0 * f64::from(called_count)))
    }

    fn total_count(self) -> Result<u64> {
        self.called_count
            .checked_add(self.missing_count)
            .ok_or_else(|| {
                GenoioError::invalid_source(
                    "<filter>",
                    "genotype count exceeds supported metadata range",
                )
            })
    }

    fn minor_allele_count(self) -> Result<Option<f64>> {
        let Some(called_alleles) = self.called_alleles()? else {
            return Ok(None);
        };
        Ok(Some(
            self.allele_count.min(called_alleles - self.allele_count),
        ))
    }

    fn missing_rate(self) -> Result<f64> {
        let total = self.total_count()?;
        if total == 0 {
            Ok(0.0)
        } else {
            Ok(self.missing_count as f64 / total as f64)
        }
    }

    fn is_polymorphic(self) -> Result<bool> {
        Ok(self.minor_allele_count()?.is_some_and(|mac| mac > 0.0))
    }

    fn mac_in_range(self, min: Option<u32>, max: Option<u32>) -> Result<bool> {
        let Some(mac) = self.minor_allele_count()? else {
            return Ok(false);
        };
        Ok(min.is_none_or(|threshold| mac >= f64::from(threshold))
            && max.is_none_or(|threshold| mac <= f64::from(threshold)))
    }

    fn maf_in_range(self, min: Option<f32>, max: Option<f32>) -> Result<bool> {
        let Some(mac) = self.minor_allele_count()? else {
            return Ok(false);
        };
        let Some(called_alleles) = self.called_alleles()? else {
            return Ok(false);
        };
        let maf = mac / called_alleles;
        Ok(min.is_none_or(|threshold| maf >= f64::from(threshold))
            && max.is_none_or(|threshold| maf <= f64::from(threshold)))
    }
}

fn dosage_counts_for_filter(values: &[f32], missing: &[bool]) -> Result<DosageFilterCounts> {
    if values.len() != missing.len() {
        return Err(GenoioError::invalid_source(
            "<filter>",
            "variant values and missing mask lengths differ",
        ));
    }

    let mut counts = DosageFilterCounts::default();
    for (value, is_missing) in values.iter().zip(missing) {
        if *is_missing {
            counts.missing_count += 1;
            continue;
        }
        if !(0.0..=2.0).contains(value) {
            return Err(GenoioError::invalid_source(
                "<filter>",
                format!("dosage statistics require values in [0, 2]; observed {value}"),
            ));
        }
        counts.allele_count += f64::from(*value);
        counts.called_count += 1;
    }
    Ok(counts)
}

fn evaluate_dosage_filter(
    values: &[f32],
    missing: &[bool],
    filter: &VariantFilter,
    variant: &VariantRecord,
    require_stats: bool,
) -> Result<(bool, Option<VariantStats>)> {
    let plan = filter.genotype_filter_plan();
    if !require_stats {
        // Matrix-only reads only need the retain/drop decision. The caller has
        // already run metadata partial evaluation, so compiled genotype plans
        // can bypass `VariantStats` construction for common dosage predicates.
        if matches!(plan, GenotypeFilterPlan::Polymorphic) {
            return Ok((is_dosage_polymorphic(values, missing)?, None));
        }
        if let Some(retain) = dosage_counts_for_filter(values, missing)?.evaluate_plan(plan)? {
            return Ok((retain, None));
        }
    }

    let stats = compute_dosage_variant_stats(values, missing)?;
    Ok((filter.evaluate(variant, Some(&stats)), Some(stats)))
}

/// Read BGEN sample and variant metadata without returning dosages.
pub fn read_bgen_metadata(bgen: &Path, sample: Option<&Path>) -> Result<MetadataOutput> {
    let mut reader = open_bgen_reader(bgen)?;
    let header = BgenHeader::read_from(&mut reader, bgen)?;
    header.validate(bgen)?;

    let samples = read_bgen_samples(&mut reader, bgen, sample, &header)?;

    reader
        .seek(SeekFrom::Start(u64::from(header.offset) + 4))
        .map_err(|source| GenoioError::Io {
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
    let mut reader = open_bgen_reader(bgen)?;
    let header = BgenHeader::read_from(&mut reader, bgen)?;
    header.validate(bgen)?;
    let all_samples = read_bgen_samples(&mut reader, bgen, sample, &header)?;
    let selection = select_samples_source_order(&all_samples, requested_samples, bgen)?;
    let mut diagnostics = selection.diagnostics.clone();
    if variant_filter.is_some_and(VariantFilter::is_always_false) {
        diagnostics.retained_variants = 0;
        return finish_dense_matrix(
            DenseMatrixParts {
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
    if matrix_only && variant_filter.is_none() {
        // Matrix-only reads do not expose variant strings or positions. Skipping
        // those bytes preserves the same matrix contract while avoiding string
        // allocation and UTF-8 validation on the hot path.
        return read_bgen_dosage_dense_matrix_only_unfiltered(
            &mut reader,
            bgen,
            &header,
            selection,
            diagnostics,
            variant_window,
        );
    }

    reader
        .seek(SeekFrom::Start(u64::from(header.offset) + 4))
        .map_err(|source| GenoioError::Io {
            path: bgen.to_path_buf(),
            source,
        })?;

    let header_variant_count = usize::try_from(header.variant_count)
        .map_err(|_| GenoioError::invalid_source(bgen, "bgen variant count is out of range"))?;
    let output_variant_capacity = variant_window.map_or(header_variant_count, |window| {
        window.len.min(header_variant_count)
    });
    let n_samples = selection.samples.len();
    let (mut values, mut missing_mask) = sample_major_buffers(n_samples, output_variant_capacity)?;
    let mut variants = Vec::with_capacity(output_variant_capacity);
    let mut decode_buffers = DosageDecodeBuffers::default();
    let mut retention = RetainedVariantState::new(variant_window);
    let mut output_variant_count = 0_usize;

    for _ in 0..header.variant_count {
        let mut variant = read_layout2_variant_identifying_data(&mut reader, bgen)?;
        let partial_decision = variant_filter.map_or(PartialFilterDecision::Accept, |filter| {
            filter.partial_decision(&variant)
        });
        match retention.metadata_decision(partial_decision, &mut diagnostics) {
            MetadataRetentionAction::Skip => {
                skip_layout2_probability_payload(&mut reader, bgen, header.flags.compression)?;
                continue;
            }
            MetadataRetentionAction::Stop => break,
            MetadataRetentionAction::Include => {
                read_layout2_probability_payload_into(
                    &mut reader,
                    bgen,
                    header.flags.compression,
                    &mut decode_buffers.payload,
                    &mut decode_buffers.compressed_payload,
                )?;
            }
            MetadataRetentionAction::DecodeGenotypes => {
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
                let (retain_variant, stats) = evaluate_dosage_filter(
                    &decode_buffers.selected_values,
                    &decode_buffers.selected_missing,
                    variant_filter.ok_or_else(|| {
                        GenoioError::internal_contract(
                            "genotype decision requires a variant filter",
                        )
                    })?,
                    &variant,
                    !matrix_only,
                )?;
                match retention.genotype_decision(retain_variant, &mut diagnostics) {
                    RetentionAction::Include => {}
                    RetentionAction::Skip => continue,
                    RetentionAction::Stop => {
                        break;
                    }
                }
                if let Some(stats) = stats {
                    attach_variant_stats(&mut variant, stats);
                }
            }
        }

        write_dosage_slot(
            BgenDosageSlotWrite {
                bgen,
                sample_count: header.sample_count,
                source_indices: &selection.source_indices,
                buffers: &mut decode_buffers,
                values: &mut values,
                missing_mask: &mut missing_mask,
                row_width: output_variant_capacity,
                variant_index: output_variant_count,
            },
            matches!(partial_decision, PartialFilterDecision::NeedGenotypes),
        )?;
        if !matrix_only {
            variants.push(variant);
        }
        output_variant_count += 1;
        if retention.window_is_satisfied() {
            break;
        }
    }

    let n_variants = output_variant_count;
    diagnostics.retained_variants = n_variants;
    shrink_sample_major_width(&mut values, n_samples, output_variant_capacity, n_variants);
    shrink_sample_major_width(
        &mut missing_mask,
        n_samples,
        output_variant_capacity,
        n_variants,
    );
    finish_dense_matrix(
        DenseMatrixParts {
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

fn open_bgen_reader(bgen: &Path) -> Result<BufReader<File>> {
    let file = File::open(bgen).map_err(|source| GenoioError::Io {
        path: bgen.to_path_buf(),
        source,
    })?;
    Ok(BufReader::with_capacity(BGEN_READER_BUFFER_SIZE, file))
}

/// Read a matrix-only BGEN block without materializing unused variant metadata.
fn read_bgen_dosage_dense_matrix_only_unfiltered(
    reader: &mut BufReader<File>,
    bgen: &Path,
    header: &BgenHeader,
    selection: DenseSampleSelection,
    mut diagnostics: genoio_core::DenseDiagnostics,
    variant_window: Option<VariantWindow>,
) -> Result<DenseGenotypeMatrix> {
    reader
        .seek(SeekFrom::Start(u64::from(header.offset) + 4))
        .map_err(|source| GenoioError::Io {
            path: bgen.to_path_buf(),
            source,
        })?;

    let header_variant_count = usize::try_from(header.variant_count)
        .map_err(|_| GenoioError::invalid_source(bgen, "bgen variant count is out of range"))?;
    let output_variant_capacity = variant_window.map_or(header_variant_count, |window| {
        window.len.min(header_variant_count)
    });
    let n_samples = selection.samples.len();
    let (mut values, mut missing_mask) = sample_major_buffers(n_samples, output_variant_capacity)?;
    let mut decode_buffers = DosageDecodeBuffers::default();
    let mut retention = RetainedVariantState::new(variant_window);
    let mut output_variant_count = 0_usize;

    for _ in 0..header.variant_count {
        if retention.window_is_satisfied() {
            break;
        }
        match retention.metadata_decision(PartialFilterDecision::Accept, &mut diagnostics) {
            MetadataRetentionAction::Include => {
                skip_layout2_variant_identifying_data(reader, bgen)?;
                read_layout2_probability_payload_into(
                    reader,
                    bgen,
                    header.flags.compression,
                    &mut decode_buffers.payload,
                    &mut decode_buffers.compressed_payload,
                )?;
                write_dosage_slot(
                    BgenDosageSlotWrite {
                        bgen,
                        sample_count: header.sample_count,
                        source_indices: &selection.source_indices,
                        buffers: &mut decode_buffers,
                        values: &mut values,
                        missing_mask: &mut missing_mask,
                        row_width: output_variant_capacity,
                        variant_index: output_variant_count,
                    },
                    false,
                )?;
                output_variant_count += 1;
            }
            MetadataRetentionAction::Skip => {
                skip_layout2_variant_identifying_data(reader, bgen)?;
                skip_layout2_probability_payload(reader, bgen, header.flags.compression)?;
            }
            MetadataRetentionAction::Stop => break,
            MetadataRetentionAction::DecodeGenotypes => {
                return Err(GenoioError::internal_contract(
                    "unfiltered bgen matrix-only path requested genotype filtering",
                ));
            }
        }
    }

    let n_variants = output_variant_count;
    diagnostics.retained_variants = n_variants;
    shrink_sample_major_width(&mut values, n_samples, output_variant_capacity, n_variants);
    shrink_sample_major_width(
        &mut missing_mask,
        n_samples,
        output_variant_capacity,
        n_variants,
    );
    finish_dense_matrix(
        DenseMatrixParts {
            n_samples,
            n_variants,
            values,
            missing_mask,
            samples: selection.samples,
            variants: Vec::new(),
            diagnostics,
        },
        true,
    )
}

/// Read retained BGEN biallelic diploid phased dosages as dense haplotype rows.
pub fn read_bgen_haplotypes_dosage_dense_windowed(
    bgen: &Path,
    sample: Option<&Path>,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    matrix_only: bool,
) -> Result<DenseGenotypeMatrix> {
    let mut reader = open_bgen_reader(bgen)?;
    let header = BgenHeader::read_from(&mut reader, bgen)?;
    header.validate(bgen)?;
    let all_samples = read_bgen_samples(&mut reader, bgen, sample, &header)?;
    let selection = select_samples_source_order(&all_samples, requested_samples, bgen)?;
    let haplotype_samples = expand_selected_samples_to_haplotypes(&selection);
    let mut diagnostics = selection.diagnostics.clone();
    if variant_filter.is_some_and(VariantFilter::is_always_false) {
        diagnostics.retained_variants = 0;
        return finish_dense_matrix(
            DenseMatrixParts {
                n_samples: haplotype_samples.len(),
                n_variants: 0,
                values: Vec::new(),
                missing_mask: Vec::new(),
                samples: haplotype_samples,
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
        return read_bgen_haplotypes_dosage_dense_indexed(context, &index_records);
    }

    reader
        .seek(SeekFrom::Start(u64::from(header.offset) + 4))
        .map_err(|source| GenoioError::Io {
            path: bgen.to_path_buf(),
            source,
        })?;

    let header_variant_count = usize::try_from(header.variant_count)
        .map_err(|_| GenoioError::invalid_source(bgen, "bgen variant count is out of range"))?;
    let output_variant_capacity = variant_window.map_or(header_variant_count, |window| {
        window.len.min(header_variant_count)
    });
    let n_haplotypes = selection.samples.len() * 2;
    let mut variants = Vec::with_capacity(output_variant_capacity);
    let mut variant_major_values = Vec::with_capacity(n_haplotypes * output_variant_capacity);
    let mut variant_major_missing = Vec::with_capacity(n_haplotypes * output_variant_capacity);
    let mut decode_buffers = HaplotypeDecodeBuffers::default();
    let mut retention = RetainedVariantState::new(variant_window);
    let mut output_variant_count = 0_usize;

    for _ in 0..header.variant_count {
        let mut variant = read_layout2_variant_identifying_data(&mut reader, bgen)?;
        let partial_decision = variant_filter.map_or(PartialFilterDecision::Accept, |filter| {
            filter.partial_decision(&variant)
        });
        match retention.metadata_decision(partial_decision, &mut diagnostics) {
            MetadataRetentionAction::Skip => {
                skip_layout2_probability_payload(&mut reader, bgen, header.flags.compression)?;
                continue;
            }
            MetadataRetentionAction::Stop => break,
            MetadataRetentionAction::Include => {
                read_layout2_probability_payload_into(
                    &mut reader,
                    bgen,
                    header.flags.compression,
                    &mut decode_buffers.payload,
                    &mut decode_buffers.compressed_payload,
                )?;
            }
            MetadataRetentionAction::DecodeGenotypes => {
                read_layout2_probability_payload_into(
                    &mut reader,
                    bgen,
                    header.flags.compression,
                    &mut decode_buffers.payload,
                    &mut decode_buffers.compressed_payload,
                )?;
                decode_buffered_haplotype_values(
                    bgen,
                    header.sample_count,
                    &selection.source_indices,
                    &mut decode_buffers,
                )?;
                let (retain_variant, stats) = evaluate_dosage_filter(
                    &decode_buffers.selected_collapsed_values,
                    &decode_buffers.selected_collapsed_missing,
                    variant_filter.ok_or_else(|| {
                        GenoioError::internal_contract(
                            "genotype decision requires a variant filter",
                        )
                    })?,
                    &variant,
                    !matrix_only,
                )?;
                match retention.genotype_decision(retain_variant, &mut diagnostics) {
                    RetentionAction::Include => {}
                    RetentionAction::Skip => continue,
                    RetentionAction::Stop => {
                        break;
                    }
                }
                if let Some(stats) = stats {
                    attach_variant_stats(&mut variant, stats);
                }
            }
        }

        if !matches!(partial_decision, PartialFilterDecision::NeedGenotypes) {
            decode_buffered_haplotype_values(
                bgen,
                header.sample_count,
                &selection.source_indices,
                &mut decode_buffers,
            )?;
        }
        if !matrix_only {
            variants.push(variant);
        }
        variant_major_values.extend_from_slice(&decode_buffers.selected_haplotype_values);
        variant_major_missing.extend_from_slice(&decode_buffers.selected_haplotype_missing);
        output_variant_count += 1;
        if retention.window_is_satisfied() {
            break;
        }
    }

    let n_samples = n_haplotypes;
    let n_variants = output_variant_count;
    diagnostics.retained_variants = n_variants;
    finish_variant_major_dense_matrix(
        VariantMajorDenseParts {
            n_samples,
            n_variants,
            variant_major_values,
            variant_major_missing,
            samples: haplotype_samples,
            variants,
            diagnostics,
        },
        matrix_only,
    )
}

fn expand_selected_samples_to_haplotypes(selection: &DenseSampleSelection) -> Vec<SampleRecord> {
    let mut haplotype_samples = Vec::with_capacity(selection.samples.len() * 2);
    for (sample, &source_index) in selection.samples.iter().zip(&selection.source_indices) {
        for haplotype_index in 0..2 {
            let mut haplotype_sample = sample.clone();
            haplotype_sample.source_sample_index = Some(source_index);
            haplotype_sample.haplotype_index = Some(haplotype_index);
            haplotype_samples.push(haplotype_sample);
        }
    }
    haplotype_samples
}

fn sample_major_buffers(n_samples: usize, n_variants: usize) -> Result<(Vec<f32>, Vec<bool>)> {
    let len = n_samples.checked_mul(n_variants).ok_or_else(|| {
        GenoioError::internal_contract("sample-major dense matrix shape is out of range")
    })?;
    Ok((vec![0.0; len], vec![false; len]))
}

struct BgenDosageSlotWrite<'a> {
    bgen: &'a Path,
    sample_count: u32,
    source_indices: &'a [usize],
    buffers: &'a mut DosageDecodeBuffers,
    values: &'a mut [f32],
    missing_mask: &'a mut [bool],
    row_width: usize,
    variant_index: usize,
}

fn write_dosage_slot(
    request: BgenDosageSlotWrite<'_>,
    already_decoded_for_filter: bool,
) -> Result<()> {
    let BgenDosageSlotWrite {
        bgen,
        sample_count,
        source_indices,
        buffers,
        values,
        missing_mask,
        row_width,
        variant_index,
    } = request;

    if !already_decoded_for_filter {
        let mut slot = SampleMajorSlotMut {
            values,
            missing_mask,
            row_width,
            variant_index,
        };
        // UKB-like unphased 8-bit records can fill the final matrix slot
        // directly. Other BGEN shapes fall back to the generic selected decode.
        if try_decode_buffered_dosage_values_into_sample_major_slot(
            bgen,
            sample_count,
            source_indices,
            buffers,
            &mut slot,
        )? {
            return Ok(());
        }
        decode_buffered_dosage_values(bgen, sample_count, source_indices, buffers)?;
    }

    write_sample_major_variant_slot(
        values,
        missing_mask,
        source_indices.len(),
        row_width,
        variant_index,
        &buffers.selected_values,
        &buffers.selected_missing,
    )
}

struct BgenIndexedReadContext<'a> {
    reader: &'a mut BufReader<File>,
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
    let n_samples = selection.samples.len();
    let (mut values, mut missing_mask) = sample_major_buffers(n_samples, output_variant_capacity)?;
    let mut variants = Vec::with_capacity(output_variant_capacity);
    let mut decode_buffers = DosageDecodeBuffers::default();
    let mut retention = RetainedVariantState::new(variant_window);
    let mut output_variant_count = 0_usize;

    for index_record in index_records {
        if retention.window_is_satisfied() {
            break;
        }
        reader
            .seek(SeekFrom::Start(index_record.file_start_position))
            .map_err(|source| GenoioError::Io {
                path: bgen.to_path_buf(),
                source,
            })?;
        let mut variant = read_layout2_variant_identifying_data(reader, bgen)?;
        let partial_decision = variant_filter.map_or(PartialFilterDecision::Accept, |filter| {
            filter.partial_decision(&variant)
        });
        match retention.metadata_decision(partial_decision, &mut diagnostics) {
            MetadataRetentionAction::Skip => {
                skip_layout2_probability_payload(reader, bgen, header.flags.compression)?;
                validate_index_record_consumed(reader, bgen, index_record)?;
                continue;
            }
            MetadataRetentionAction::Stop => {
                skip_layout2_probability_payload(reader, bgen, header.flags.compression)?;
                validate_index_record_consumed(reader, bgen, index_record)?;
                break;
            }
            MetadataRetentionAction::Include => {
                read_layout2_probability_payload_into(
                    reader,
                    bgen,
                    header.flags.compression,
                    &mut decode_buffers.payload,
                    &mut decode_buffers.compressed_payload,
                )?;
            }
            MetadataRetentionAction::DecodeGenotypes => {
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
                let (retain_variant, stats) = evaluate_dosage_filter(
                    &decode_buffers.selected_values,
                    &decode_buffers.selected_missing,
                    variant_filter.ok_or_else(|| {
                        GenoioError::internal_contract(
                            "genotype decision requires a variant filter",
                        )
                    })?,
                    &variant,
                    !matrix_only,
                )?;
                match retention.genotype_decision(retain_variant, &mut diagnostics) {
                    RetentionAction::Include => {}
                    RetentionAction::Skip => {
                        validate_index_record_consumed(reader, bgen, index_record)?;
                        continue;
                    }
                    RetentionAction::Stop => {
                        validate_index_record_consumed(reader, bgen, index_record)?;
                        break;
                    }
                }
                if let Some(stats) = stats {
                    attach_variant_stats(&mut variant, stats);
                }
            }
        }

        write_dosage_slot(
            BgenDosageSlotWrite {
                bgen,
                sample_count: header.sample_count,
                source_indices: &selection.source_indices,
                buffers: &mut decode_buffers,
                values: &mut values,
                missing_mask: &mut missing_mask,
                row_width: output_variant_capacity,
                variant_index: output_variant_count,
            },
            matches!(partial_decision, PartialFilterDecision::NeedGenotypes),
        )?;
        validate_index_record_consumed(reader, bgen, index_record)?;
        if !matrix_only {
            variants.push(variant);
        }
        output_variant_count += 1;
    }

    let n_variants = output_variant_count;
    diagnostics.retained_variants = n_variants;
    shrink_sample_major_width(&mut values, n_samples, output_variant_capacity, n_variants);
    shrink_sample_major_width(
        &mut missing_mask,
        n_samples,
        output_variant_capacity,
        n_variants,
    );
    finish_dense_matrix(
        DenseMatrixParts {
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

fn read_bgen_haplotypes_dosage_dense_indexed(
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
    let haplotype_samples = expand_selected_samples_to_haplotypes(&selection);
    let output_variant_capacity = variant_window.map_or(index_records.len(), |window| {
        window.len.min(index_records.len())
    });
    let n_haplotypes = selection.samples.len() * 2;
    let mut variants = Vec::with_capacity(output_variant_capacity);
    let mut variant_major_values = Vec::with_capacity(n_haplotypes * output_variant_capacity);
    let mut variant_major_missing = Vec::with_capacity(n_haplotypes * output_variant_capacity);
    let mut decode_buffers = HaplotypeDecodeBuffers::default();
    let mut retention = RetainedVariantState::new(variant_window);
    let mut output_variant_count = 0_usize;

    for index_record in index_records {
        if retention.window_is_satisfied() {
            break;
        }
        reader
            .seek(SeekFrom::Start(index_record.file_start_position))
            .map_err(|source| GenoioError::Io {
                path: bgen.to_path_buf(),
                source,
            })?;
        let mut variant = read_layout2_variant_identifying_data(reader, bgen)?;
        let partial_decision = variant_filter.map_or(PartialFilterDecision::Accept, |filter| {
            filter.partial_decision(&variant)
        });
        match retention.metadata_decision(partial_decision, &mut diagnostics) {
            MetadataRetentionAction::Skip => {
                skip_layout2_probability_payload(reader, bgen, header.flags.compression)?;
                validate_index_record_consumed(reader, bgen, index_record)?;
                continue;
            }
            MetadataRetentionAction::Stop => {
                skip_layout2_probability_payload(reader, bgen, header.flags.compression)?;
                validate_index_record_consumed(reader, bgen, index_record)?;
                break;
            }
            MetadataRetentionAction::Include => {
                read_layout2_probability_payload_into(
                    reader,
                    bgen,
                    header.flags.compression,
                    &mut decode_buffers.payload,
                    &mut decode_buffers.compressed_payload,
                )?;
            }
            MetadataRetentionAction::DecodeGenotypes => {
                read_layout2_probability_payload_into(
                    reader,
                    bgen,
                    header.flags.compression,
                    &mut decode_buffers.payload,
                    &mut decode_buffers.compressed_payload,
                )?;
                decode_buffered_haplotype_values(
                    bgen,
                    header.sample_count,
                    &selection.source_indices,
                    &mut decode_buffers,
                )?;
                let (retain_variant, stats) = evaluate_dosage_filter(
                    &decode_buffers.selected_collapsed_values,
                    &decode_buffers.selected_collapsed_missing,
                    variant_filter.ok_or_else(|| {
                        GenoioError::internal_contract(
                            "genotype decision requires a variant filter",
                        )
                    })?,
                    &variant,
                    !matrix_only,
                )?;
                match retention.genotype_decision(retain_variant, &mut diagnostics) {
                    RetentionAction::Include => {}
                    RetentionAction::Skip => {
                        validate_index_record_consumed(reader, bgen, index_record)?;
                        continue;
                    }
                    RetentionAction::Stop => {
                        validate_index_record_consumed(reader, bgen, index_record)?;
                        break;
                    }
                }
                if let Some(stats) = stats {
                    attach_variant_stats(&mut variant, stats);
                }
            }
        }

        if !matches!(partial_decision, PartialFilterDecision::NeedGenotypes) {
            decode_buffered_haplotype_values(
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
        variant_major_values.extend_from_slice(&decode_buffers.selected_haplotype_values);
        variant_major_missing.extend_from_slice(&decode_buffers.selected_haplotype_missing);
        output_variant_count += 1;
    }

    let n_samples = n_haplotypes;
    let n_variants = output_variant_count;
    diagnostics.retained_variants = n_variants;
    finish_variant_major_dense_matrix(
        VariantMajorDenseParts {
            n_samples,
            n_variants,
            variant_major_values,
            variant_major_missing,
            samples: haplotype_samples,
            variants,
            diagnostics,
        },
        matrix_only,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use genoio_core::GenotypeFilterConjunction;
    use serde_json::json;

    fn test_variant() -> VariantRecord {
        VariantRecord {
            chrom: "1".to_string(),
            pos: 10,
            id: "rs1".to_string(),
            a0: "A".to_string(),
            a1: "G".to_string(),
            ref_allele: Some("A".to_string()),
            alt_allele: Some("G".to_string()),
            source_a0: "A".to_string(),
            source_a1: "G".to_string(),
            flipped: false,
            qual: None,
            af: None,
            maf: None,
            mac: None,
            missing_rate: None,
            n_called: None,
        }
    }

    fn genotype_filter(name: &str, params: serde_json::Value) -> VariantFilter {
        VariantFilter::from_json_value(json!({
            "op": "predicate",
            "name": name,
            "params": params,
        }))
        .unwrap()
    }

    fn dosage_fixture() -> ([f32; 4], [bool; 4]) {
        ([0.0, 1.0, 2.0, 2.0], [false, false, false, true])
    }

    #[test]
    fn dosage_filter_plan_evaluates_mac_maf_and_missing_rate() {
        let (values, missing) = dosage_fixture();
        let counts = dosage_counts_for_filter(&values, &missing).unwrap();

        assert_eq!(
            counts
                .evaluate_plan(GenotypeFilterPlan::MacRange {
                    min: Some(3),
                    max: Some(3),
                })
                .unwrap(),
            Some(true)
        );
        assert_eq!(
            counts
                .evaluate_plan(GenotypeFilterPlan::MafRange {
                    min: Some(0.49),
                    max: Some(0.51),
                })
                .unwrap(),
            Some(true)
        );
        assert_eq!(
            counts
                .evaluate_plan(GenotypeFilterPlan::MissingRateMax { max: 0.20 })
                .unwrap(),
            Some(false)
        );
    }

    #[test]
    fn dosage_filter_plan_evaluates_conjunctions() {
        let (values, missing) = dosage_fixture();
        let counts = dosage_counts_for_filter(&values, &missing).unwrap();

        let passing = GenotypeFilterConjunction {
            polymorphic: true,
            mac_min: Some(2),
            mac_max: Some(4),
            maf_min: Some(0.4),
            maf_max: Some(0.6),
            missing_rate_max: Some(0.3),
        };
        let failing = GenotypeFilterConjunction {
            missing_rate_max: Some(0.2),
            ..passing
        };

        assert_eq!(
            counts
                .evaluate_plan(GenotypeFilterPlan::Conjunction(passing))
                .unwrap(),
            Some(true)
        );
        assert_eq!(
            counts
                .evaluate_plan(GenotypeFilterPlan::Conjunction(failing))
                .unwrap(),
            Some(false)
        );
    }

    #[test]
    fn matrix_only_dosage_filter_returns_decision_without_stats() {
        let (values, missing) = dosage_fixture();
        let filter = genotype_filter("maf", json!({ "min": 0.49, "max": 0.51 }));

        let (retain, stats) =
            evaluate_dosage_filter(&values, &missing, &filter, &test_variant(), false).unwrap();

        assert!(retain);
        assert_eq!(stats, None);
    }

    #[test]
    fn metadata_output_dosage_filter_keeps_stats() {
        let (values, missing) = dosage_fixture();
        let filter = genotype_filter("maf", json!({ "min": 0.49, "max": 0.51 }));

        let (retain, stats) =
            evaluate_dosage_filter(&values, &missing, &filter, &test_variant(), true).unwrap();

        assert!(retain);
        assert_eq!(stats.unwrap().missing_rate, 0.25);
    }

    #[test]
    fn dosage_filter_counts_match_variant_stats_thresholds() {
        let (values, missing) = dosage_fixture();
        let counts = dosage_counts_for_filter(&values, &missing).unwrap();
        let stats = compute_dosage_variant_stats(&values, &missing).unwrap();

        assert_eq!(counts.is_polymorphic().unwrap(), stats.polymorphic);
        assert_eq!(counts.missing_rate().unwrap(), stats.missing_rate);
        assert_eq!(counts.minor_allele_count().unwrap(), stats.mac);
    }
}
