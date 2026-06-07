// pattern: Imperative Shell
//! BGEN reader orchestration and matrix assembly.

use std::fs::File;
use std::io::{Seek, SeekFrom};
use std::path::Path;

use crate::matrix::{
    finish_dense_matrix, finish_variant_major_dense_matrix, DenseMatrixParts,
    VariantMajorDenseParts,
};
use crate::retention::{MetadataRetentionAction, RetainedVariantState, RetentionAction};
use crate::Result;
use genoio_core::{
    attach_variant_stats, compute_dosage_variant_stats, select_samples_source_order,
    DenseGenotypeMatrix, DenseSampleSelection, GenoioError, MetadataOutput, PartialFilterDecision,
    SampleRecord, SourceCapabilities, VariantFilter, VariantWindow,
};

mod decode;
mod header;
mod index;
mod io;

use decode::{
    decode_buffered_dosage_values, decode_buffered_haplotype_values,
    read_layout2_probability_payload_into, skip_layout2_probability_payload, DosageDecodeBuffers,
    HaplotypeDecodeBuffers,
};
use header::{
    read_bgen_samples, read_layout2_variant_identifying_data, read_layout2_variant_metadata,
    BgenHeader,
};
use index::{indexed_region_records, validate_index_record_consumed, BgenIndexRecord};

/// Read BGEN sample and variant metadata without returning dosages.
pub fn read_bgen_metadata(bgen: &Path, sample: Option<&Path>) -> Result<MetadataOutput> {
    let mut reader = File::open(bgen).map_err(|source| GenoioError::Io {
        path: bgen.to_path_buf(),
        source,
    })?;
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
    let mut reader = File::open(bgen).map_err(|source| GenoioError::Io {
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
    let mut variants = Vec::with_capacity(output_variant_capacity);
    let mut variant_major_values =
        Vec::with_capacity(selection.samples.len() * output_variant_capacity);
    let mut variant_major_missing =
        Vec::with_capacity(selection.samples.len() * output_variant_capacity);
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
                let stats = compute_dosage_variant_stats(
                    &decode_buffers.selected_values,
                    &decode_buffers.selected_missing,
                )?;
                match retention.genotype_decision(
                    variant_filter.is_none_or(|filter| filter.evaluate(&variant, Some(&stats))),
                    &mut diagnostics,
                ) {
                    RetentionAction::Include => {}
                    RetentionAction::Skip => continue,
                    RetentionAction::Stop => {
                        break;
                    }
                }
                attach_variant_stats(&mut variant, stats);
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
        if retention.window_is_satisfied() {
            break;
        }
    }

    let n_samples = selection.samples.len();
    let n_variants = output_variant_count;
    diagnostics.retained_variants = n_variants;
    finish_variant_major_dense_matrix(
        VariantMajorDenseParts {
            n_samples,
            n_variants,
            variant_major_values,
            variant_major_missing,
            samples: selection.samples,
            variants,
            diagnostics,
        },
        matrix_only,
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
    let mut reader = File::open(bgen).map_err(|source| GenoioError::Io {
        path: bgen.to_path_buf(),
        source,
    })?;
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
                let stats = compute_dosage_variant_stats(
                    &decode_buffers.selected_collapsed_values,
                    &decode_buffers.selected_collapsed_missing,
                )?;
                match retention.genotype_decision(
                    variant_filter.is_none_or(|filter| filter.evaluate(&variant, Some(&stats))),
                    &mut diagnostics,
                ) {
                    RetentionAction::Include => {}
                    RetentionAction::Skip => continue,
                    RetentionAction::Stop => {
                        break;
                    }
                }
                attach_variant_stats(&mut variant, stats);
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
                let stats = compute_dosage_variant_stats(
                    &decode_buffers.selected_values,
                    &decode_buffers.selected_missing,
                )?;
                match retention.genotype_decision(
                    variant_filter.is_none_or(|filter| filter.evaluate(&variant, Some(&stats))),
                    &mut diagnostics,
                ) {
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
                attach_variant_stats(&mut variant, stats);
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
    finish_variant_major_dense_matrix(
        VariantMajorDenseParts {
            n_samples,
            n_variants,
            variant_major_values,
            variant_major_missing,
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
                let stats = compute_dosage_variant_stats(
                    &decode_buffers.selected_collapsed_values,
                    &decode_buffers.selected_collapsed_missing,
                )?;
                match retention.genotype_decision(
                    variant_filter.is_none_or(|filter| filter.evaluate(&variant, Some(&stats))),
                    &mut diagnostics,
                ) {
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
                attach_variant_stats(&mut variant, stats);
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
