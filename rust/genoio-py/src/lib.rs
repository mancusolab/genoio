// pattern: Imperative Shell

use std::path::PathBuf;
use std::slice;

use pyo3::prelude::*;
use pyo3::types::{
    PyBool, PyByteArray, PyDict, PyFloat, PyInt, PyList, PyModule, PyString, PyTuple,
};

#[pyfunction]
/// Return the Rust IO backend name for Python diagnostics.
fn backend_name() -> &'static str {
    genoio_io::backend_name()
}

#[pyfunction]
/// Read source metadata and convert it into Python dictionaries.
fn read_metadata(
    py: Python<'_>,
    format: &str,
    members: &Bound<'_, PyDict>,
) -> PyResult<Py<PyDict>> {
    let output = match format {
        "vcf" | "bcf" => {
            let key = if format == "vcf" { "vcf" } else { "bcf" };
            let path = member_path(members, key)?;
            genoio_io::read_vcf_metadata(&path)
        }
        "plink1" => {
            let bed = member_path(members, "bed")?;
            let bim = member_path(members, "bim")?;
            let fam = member_path(members, "fam")?;
            genoio_io::read_plink1_metadata(&bed, &bim, &fam)
        }
        "plink2" => {
            let pgen = member_path(members, "pgen")?;
            let pvar = member_path(members, "pvar")?;
            let psam = member_path(members, "psam")?;
            genoio_io::read_plink2_metadata(&pgen, &pvar, &psam)
        }
        "bgen" => {
            let bgen = member_path(members, "bgen")?;
            let sample = optional_member_path(members, "sample")?;
            genoio_io::read_bgen_metadata(&bgen, sample.as_deref())
        }
        other => {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "unsupported metadata format: {other}"
            )));
        }
    }
    .map_err(|error| pyo3::exceptions::PyValueError::new_err(error.to_string()))?;

    metadata_to_py(py, output)
}

#[pyfunction]
/// Read dense genotypes from a resolved source.
fn read_dense(
    py: Python<'_>,
    format: &str,
    members: &Bound<'_, PyDict>,
    options: &Bound<'_, PyDict>,
) -> PyResult<Py<PyDict>> {
    let read_options = read_options(options)?;
    let output = match format {
        "vcf" | "bcf" => {
            let key = if format == "vcf" { "vcf" } else { "bcf" };
            let path = member_path(members, key)?;
            match read_options.dosage {
                DosageSource::Hardcall => genoio_io::read_vcf_dense_windowed(
                    &path,
                    read_options.requested_samples.as_deref(),
                    read_options.variant_filter.as_ref(),
                    read_options.variant_window,
                ),
                DosageSource::Dosage => genoio_io::read_vcf_dosage_dense_windowed(
                    &path,
                    read_options.requested_samples.as_deref(),
                    read_options.variant_filter.as_ref(),
                    read_options.variant_window,
                ),
            }
        }
        "plink1" => {
            if read_options.dosage != DosageSource::Hardcall {
                return Err(pyo3::exceptions::PyValueError::new_err(
                    "plink1 does not support dosage-backed genotype reads",
                ));
            }
            let bed = member_path(members, "bed")?;
            let bim = member_path(members, "bim")?;
            let fam = member_path(members, "fam")?;
            genoio_io::read_plink1_dense_windowed(
                &bed,
                &bim,
                &fam,
                read_options.requested_samples.as_deref(),
                read_options.variant_filter.as_ref(),
                read_options.variant_window,
            )
        }
        "plink2" => {
            let pgen = member_path(members, "pgen")?;
            let pvar = member_path(members, "pvar")?;
            let psam = member_path(members, "psam")?;
            match read_options.dosage {
                DosageSource::Hardcall => genoio_io::read_plink2_dense_windowed(
                    &pgen,
                    &pvar,
                    &psam,
                    read_options.requested_samples.as_deref(),
                    read_options.variant_filter.as_ref(),
                    read_options.variant_window,
                    read_options.matrix_only,
                ),
                DosageSource::Dosage => genoio_io::read_plink2_dosage_dense_windowed(
                    &pgen,
                    &pvar,
                    &psam,
                    read_options.requested_samples.as_deref(),
                    read_options.variant_filter.as_ref(),
                    read_options.variant_window,
                ),
            }
        }
        "bgen" => {
            let bgen = member_path(members, "bgen")?;
            let sample = optional_member_path(members, "sample")?;
            match read_options.dosage {
                DosageSource::Hardcall => {
                    return Err(pyo3::exceptions::PyValueError::new_err(
                        "bgen hardcall genotype reads are not implemented",
                    ));
                }
                DosageSource::Dosage => genoio_io::read_bgen_dosage_dense_windowed(
                    &bgen,
                    sample.as_deref(),
                    read_options.requested_samples.as_deref(),
                    read_options.variant_filter.as_ref(),
                    read_options.variant_window,
                ),
            }
        }
        other => {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "unsupported dense format: {other}"
            )));
        }
    }
    .map_err(|error| pyo3::exceptions::PyValueError::new_err(error.to_string()))?;

    dense_to_py(
        py,
        output,
        read_options.return_samples,
        read_options.return_variants,
    )
}

#[pyfunction]
/// Read sparse genotypes from a resolved source.
fn read_sparse(
    py: Python<'_>,
    format: &str,
    members: &Bound<'_, PyDict>,
    options: &Bound<'_, PyDict>,
) -> PyResult<Py<PyDict>> {
    let read_options = read_options(options)?;
    if read_options.dosage != DosageSource::Hardcall {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "sparse dosage-backed genotype reads are not implemented",
        ));
    }
    let output = match format {
        "vcf" | "bcf" => {
            let key = if format == "vcf" { "vcf" } else { "bcf" };
            let path = member_path(members, key)?;
            genoio_io::read_vcf_sparse_windowed(
                &path,
                read_options.requested_samples.as_deref(),
                read_options.variant_filter.as_ref(),
                read_options.variant_window,
            )
        }
        "plink1" => {
            let bed = member_path(members, "bed")?;
            let bim = member_path(members, "bim")?;
            let fam = member_path(members, "fam")?;
            genoio_io::read_plink1_sparse_windowed(
                &bed,
                &bim,
                &fam,
                read_options.requested_samples.as_deref(),
                read_options.variant_filter.as_ref(),
                read_options.variant_window,
            )
        }
        "plink2" => {
            let pgen = member_path(members, "pgen")?;
            let pvar = member_path(members, "pvar")?;
            let psam = member_path(members, "psam")?;
            genoio_io::read_plink2_sparse_windowed(
                &pgen,
                &pvar,
                &psam,
                read_options.requested_samples.as_deref(),
                read_options.variant_filter.as_ref(),
                read_options.variant_window,
            )
        }
        "bgen" => {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "bgen sparse genotype reads are not implemented",
            ));
        }
        other => {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "unsupported sparse format: {other}"
            )));
        }
    }
    .map_err(|error| pyo3::exceptions::PyValueError::new_err(error.to_string()))?;

    sparse_to_py(
        py,
        output,
        read_options.return_samples,
        read_options.return_variants,
    )
}

#[pyfunction]
/// Read phased VCF genotypes as dense haplotype rows.
fn read_haplotypes_dense(
    py: Python<'_>,
    format: &str,
    members: &Bound<'_, PyDict>,
    options: &Bound<'_, PyDict>,
) -> PyResult<Py<PyDict>> {
    let read_options = read_options(options)?;
    let output = match format {
        "vcf" | "bcf" => {
            let key = if format == "vcf" { "vcf" } else { "bcf" };
            let path = member_path(members, key)?;
            genoio_io::read_vcf_haplotypes_dense_windowed(
                &path,
                read_options.requested_samples.as_deref(),
                read_options.variant_filter.as_ref(),
                read_options.variant_window,
            )
        }
        "bgen" => {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "bgen does not support haplo reads",
            ));
        }
        other => {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "unsupported haplotype format: {other}"
            )));
        }
    }
    .map_err(|error| pyo3::exceptions::PyValueError::new_err(error.to_string()))?;

    dense_to_py(
        py,
        output,
        read_options.return_samples,
        read_options.return_variants,
    )
}

#[pyfunction]
/// Read phased VCF genotypes as sparse haplotype rows.
fn read_haplotypes_sparse(
    py: Python<'_>,
    format: &str,
    members: &Bound<'_, PyDict>,
    options: &Bound<'_, PyDict>,
) -> PyResult<Py<PyDict>> {
    let read_options = read_options(options)?;
    let output = match format {
        "vcf" | "bcf" => {
            let key = if format == "vcf" { "vcf" } else { "bcf" };
            let path = member_path(members, key)?;
            genoio_io::read_vcf_haplotypes_sparse_windowed(
                &path,
                read_options.requested_samples.as_deref(),
                read_options.variant_filter.as_ref(),
                read_options.variant_window,
            )
        }
        "bgen" => {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "bgen does not support haplo reads",
            ));
        }
        other => {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "unsupported haplotype format: {other}"
            )));
        }
    }
    .map_err(|error| pyo3::exceptions::PyValueError::new_err(error.to_string()))?;

    sparse_to_py(
        py,
        output,
        read_options.return_samples,
        read_options.return_variants,
    )
}

#[pymodule]
/// PyO3 extension module used by the pure-Python public API.
fn _rust(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(backend_name, module)?)?;
    module.add_function(wrap_pyfunction!(read_metadata, module)?)?;
    module.add_function(wrap_pyfunction!(read_dense, module)?)?;
    module.add_function(wrap_pyfunction!(read_sparse, module)?)?;
    module.add_function(wrap_pyfunction!(read_haplotypes_dense, module)?)?;
    module.add_function(wrap_pyfunction!(read_haplotypes_sparse, module)?)?;
    Ok(())
}

fn member_path(members: &Bound<'_, PyDict>, key: &str) -> PyResult<PathBuf> {
    let value = members.get_item(key)?.ok_or_else(|| {
        pyo3::exceptions::PyKeyError::new_err(format!("missing source member: {key}"))
    })?;
    Ok(PathBuf::from(value.extract::<String>()?))
}

fn optional_member_path(members: &Bound<'_, PyDict>, key: &str) -> PyResult<Option<PathBuf>> {
    let Some(value) = members.get_item(key)? else {
        return Ok(None);
    };
    Ok(Some(PathBuf::from(value.extract::<String>()?)))
}

struct ReadOptions {
    requested_samples: Option<Vec<String>>,
    variant_filter: Option<genoio_core::VariantFilter>,
    variant_window: Option<genoio_core::VariantWindow>,
    dosage: DosageSource,
    return_samples: bool,
    return_variants: bool,
    matrix_only: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum DosageSource {
    Hardcall,
    Dosage,
}

fn read_options(options: &Bound<'_, PyDict>) -> PyResult<ReadOptions> {
    Ok(ReadOptions {
        requested_samples: samples_option(options)?,
        variant_filter: variants_option(options)?,
        variant_window: variant_window_option(options)?,
        dosage: dosage_option(options)?,
        return_samples: bool_option(options, "return_samples")?,
        return_variants: bool_option(options, "return_variants")?,
        matrix_only: required_bool_option(options, "matrix_only")?,
    })
}

fn samples_option(options: &Bound<'_, PyDict>) -> PyResult<Option<Vec<String>>> {
    let Some(value) = options.get_item("samples")? else {
        return Ok(None);
    };
    if value.is_none() {
        return Ok(None);
    }
    value.extract::<Vec<String>>().map(Some)
}

fn variants_option(options: &Bound<'_, PyDict>) -> PyResult<Option<genoio_core::VariantFilter>> {
    let Some(value) = options.get_item("variants")? else {
        return Ok(None);
    };
    if value.is_none() {
        return Ok(None);
    }
    // Python validates and serializes the public FilterExpr tree; Rust parses
    // the JSON-compatible IR again so malformed direct extension calls fail at
    // the same boundary as normal reads.
    let json = py_to_json_value(&value)?;
    genoio_core::VariantFilter::from_json_value(json)
        .map(Some)
        .map_err(|error| pyo3::exceptions::PyValueError::new_err(error.to_string()))
}

fn variant_window_option(
    options: &Bound<'_, PyDict>,
) -> PyResult<Option<genoio_core::VariantWindow>> {
    let Some(value) = options.get_item("variant_window")? else {
        return Ok(None);
    };
    if value.is_none() {
        return Ok(None);
    }
    let dict = value.downcast::<PyDict>()?;
    let start = dict
        .get_item("start")?
        .ok_or_else(|| pyo3::exceptions::PyKeyError::new_err("missing variant_window.start"))?
        .extract::<usize>()?;
    let len = dict
        .get_item("len")?
        .ok_or_else(|| pyo3::exceptions::PyKeyError::new_err("missing variant_window.len"))?
        .extract::<usize>()?;
    Ok(Some(genoio_core::VariantWindow { start, len }))
}

fn dosage_option(options: &Bound<'_, PyDict>) -> PyResult<DosageSource> {
    let value = options
        .get_item("dosage")?
        .ok_or_else(|| pyo3::exceptions::PyKeyError::new_err("missing option: dosage"))?;
    match value.extract::<String>()?.as_str() {
        "hardcall" => Ok(DosageSource::Hardcall),
        "dosage" => Ok(DosageSource::Dosage),
        other => Err(pyo3::exceptions::PyValueError::new_err(format!(
            "unsupported dosage source: {other}"
        ))),
    }
}

fn bool_option(options: &Bound<'_, PyDict>, key: &str) -> PyResult<bool> {
    let Some(value) = options.get_item(key)? else {
        return Ok(false);
    };
    if value.is_none() {
        return Ok(false);
    }
    value.extract::<bool>()
}

fn required_bool_option(options: &Bound<'_, PyDict>, key: &str) -> PyResult<bool> {
    let value = options
        .get_item(key)?
        .ok_or_else(|| pyo3::exceptions::PyKeyError::new_err(format!("missing option: {key}")))?;
    if value.is_none() {
        return Err(pyo3::exceptions::PyTypeError::new_err(format!(
            "{key} must be a bool"
        )));
    }
    value.extract::<bool>()
}

fn metadata_to_py(py: Python<'_>, output: genoio_core::MetadataOutput) -> PyResult<Py<PyDict>> {
    let genoio_core::MetadataOutput {
        samples,
        variants,
        capabilities: source_capabilities,
    } = output;
    let dict = PyDict::new(py);

    let capabilities = PyDict::new(py);
    capabilities.set_item("supports_geno", source_capabilities.supports_geno)?;
    capabilities.set_item("supports_haplo", source_capabilities.supports_haplo)?;
    capabilities.set_item("phased", source_capabilities.phased)?;

    dict.set_item("samples", sample_records_to_py(py, samples)?)?;
    dict.set_item("variants", variant_records_to_py(py, variants)?)?;
    dict.set_item("capabilities", capabilities)?;
    Ok(dict.unbind())
}

fn dense_to_py(
    py: Python<'_>,
    output: genoio_core::DenseGenotypeMatrix,
    return_samples: bool,
    return_variants: bool,
) -> PyResult<Py<PyDict>> {
    let dict = PyDict::new(py);
    dict.set_item("values", f32_vec_to_numpy(py, output.values)?)?;
    dict.set_item("shape", (output.n_samples, output.n_variants))?;
    dict.set_item("missing_mask", bool_vec_to_numpy(py, output.missing_mask)?)?;
    if return_samples {
        dict.set_item("samples", sample_records_to_py(py, output.samples)?)?;
    }
    if return_variants {
        dict.set_item("variants", variant_records_to_py(py, output.variants)?)?;
    }

    let diagnostics = PyDict::new(py);
    diagnostics.set_item("requested_samples", output.diagnostics.requested_samples)?;
    diagnostics.set_item("retained_samples", output.diagnostics.retained_samples)?;
    diagnostics.set_item("missing_samples", output.diagnostics.missing_samples)?;
    diagnostics.set_item("candidate_variants", output.diagnostics.candidate_variants)?;
    diagnostics.set_item("retained_variants", output.diagnostics.retained_variants)?;
    diagnostics.set_item(
        "dropped_metadata_variants",
        output.diagnostics.dropped_metadata_variants,
    )?;
    diagnostics.set_item(
        "dropped_genotype_variants",
        output.diagnostics.dropped_genotype_variants,
    )?;
    dict.set_item("diagnostics", diagnostics)?;

    Ok(dict.unbind())
}

fn sparse_to_py(
    py: Python<'_>,
    output: genoio_core::SparseGenotypeMatrix,
    return_samples: bool,
    return_variants: bool,
) -> PyResult<Py<PyDict>> {
    let dict = PyDict::new(py);
    dict.set_item("indptr", usize_vec_to_numpy_i64(py, output.indptr)?)?;
    dict.set_item("indices", usize_vec_to_numpy_i64(py, output.indices)?)?;
    dict.set_item("data", f32_vec_to_numpy(py, output.data)?)?;
    dict.set_item("shape", (output.n_rows, output.n_cols))?;
    if return_samples {
        dict.set_item("samples", sample_records_to_py(py, output.samples)?)?;
    }
    if return_variants {
        dict.set_item("variants", variant_records_to_py(py, output.variants)?)?;
    }

    let diagnostics = PyDict::new(py);
    diagnostics.set_item("requested_samples", output.diagnostics.requested_samples)?;
    diagnostics.set_item("retained_samples", output.diagnostics.retained_samples)?;
    diagnostics.set_item("missing_samples", output.diagnostics.missing_samples)?;
    diagnostics.set_item("candidate_variants", output.diagnostics.candidate_variants)?;
    diagnostics.set_item("retained_variants", output.diagnostics.retained_variants)?;
    diagnostics.set_item(
        "dropped_metadata_variants",
        output.diagnostics.dropped_metadata_variants,
    )?;
    diagnostics.set_item(
        "dropped_genotype_variants",
        output.diagnostics.dropped_genotype_variants,
    )?;
    dict.set_item("diagnostics", diagnostics)?;

    Ok(dict.unbind())
}

fn f32_vec_to_numpy(py: Python<'_>, values: Vec<f32>) -> PyResult<Bound<'_, PyAny>> {
    // SAFETY: f32 has a stable byte representation for NumPy's float32 view.
    // PyByteArray owns a copy of the bytes before `values` is dropped, so the
    // returned array does not borrow Rust memory.
    let bytes = unsafe {
        slice::from_raw_parts(
            values.as_ptr().cast::<u8>(),
            values.len() * std::mem::size_of::<f32>(),
        )
    };
    let buffer = PyByteArray::new(py, bytes);
    PyModule::import(py, "numpy")?.call_method1("frombuffer", (buffer, "float32"))
}

fn bool_vec_to_numpy(py: Python<'_>, values: Vec<bool>) -> PyResult<Bound<'_, PyAny>> {
    let bytes = values.into_iter().map(u8::from).collect::<Vec<u8>>();
    let buffer = PyByteArray::new(py, &bytes);
    PyModule::import(py, "numpy")?.call_method1("frombuffer", (buffer, "bool"))
}

fn usize_vec_to_numpy_i64(py: Python<'_>, values: Vec<usize>) -> PyResult<Bound<'_, PyAny>> {
    let values = values
        .into_iter()
        .map(|value| {
            i64::try_from(value).map_err(|_| {
                pyo3::exceptions::PyOverflowError::new_err(
                    "array index exceeds supported NumPy int64 range",
                )
            })
        })
        .collect::<PyResult<Vec<i64>>>()?;
    // SAFETY: values has been converted to i64 and PyByteArray owns a copy of
    // the contiguous bytes before the local Vec is dropped.
    let bytes = unsafe {
        slice::from_raw_parts(
            values.as_ptr().cast::<u8>(),
            values.len() * std::mem::size_of::<i64>(),
        )
    };
    let buffer = PyByteArray::new(py, bytes);
    PyModule::import(py, "numpy")?.call_method1("frombuffer", (buffer, "int64"))
}

fn sample_records_to_py(
    py: Python<'_>,
    samples: Vec<genoio_core::SampleRecord>,
) -> PyResult<Bound<'_, PyDict>> {
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
        source_sample_indices.push(sample.source_sample_index);
        haplotype_indices.push(sample.haplotype_index);
    }

    let columns = PyDict::new(py);
    columns.set_item("fid", fids)?;
    columns.set_item("iid", iids)?;
    columns.set_item("father", fathers)?;
    columns.set_item("mother", mothers)?;
    columns.set_item("sex", sexes)?;
    columns.set_item("phenotype", phenotypes)?;
    columns.set_item("source_sample_index", source_sample_indices)?;
    columns.set_item("haplotype_index", haplotype_indices)?;
    Ok(columns)
}

fn variant_records_to_py(
    py: Python<'_>,
    variants: Vec<genoio_core::VariantRecord>,
) -> PyResult<Bound<'_, PyDict>> {
    let mut chroms = Vec::with_capacity(variants.len());
    let mut positions = Vec::with_capacity(variants.len());
    let mut ids = Vec::with_capacity(variants.len());
    let mut a0s = Vec::with_capacity(variants.len());
    let mut a1s = Vec::with_capacity(variants.len());
    let mut ref_alleles = Vec::with_capacity(variants.len());
    let mut alt_alleles = Vec::with_capacity(variants.len());
    let mut source_a0s = Vec::with_capacity(variants.len());
    let mut source_a1s = Vec::with_capacity(variants.len());
    let mut flipped = Vec::with_capacity(variants.len());
    let mut quals = Vec::with_capacity(variants.len());
    let mut afs = Vec::with_capacity(variants.len());
    let mut mafs = Vec::with_capacity(variants.len());
    let mut macs = Vec::with_capacity(variants.len());
    let mut missing_rates = Vec::with_capacity(variants.len());
    let mut n_called = Vec::with_capacity(variants.len());

    for variant in variants {
        chroms.push(variant.chrom);
        positions.push(variant.pos);
        ids.push(variant.id);
        a0s.push(variant.a0);
        a1s.push(variant.a1);
        ref_alleles.push(variant.ref_allele);
        alt_alleles.push(variant.alt_allele);
        source_a0s.push(variant.source_a0);
        source_a1s.push(variant.source_a1);
        flipped.push(variant.flipped);
        quals.push(variant.qual);
        afs.push(variant.af);
        mafs.push(variant.maf);
        macs.push(variant.mac);
        missing_rates.push(variant.missing_rate);
        n_called.push(variant.n_called);
    }

    let columns = PyDict::new(py);
    columns.set_item("chrom", chroms)?;
    columns.set_item("pos", positions)?;
    columns.set_item("id", ids)?;
    columns.set_item("a0", a0s)?;
    columns.set_item("a1", a1s)?;
    columns.set_item("ref_allele", ref_alleles)?;
    columns.set_item("alt_allele", alt_alleles)?;
    columns.set_item("source_a0", source_a0s)?;
    columns.set_item("source_a1", source_a1s)?;
    columns.set_item("flipped", flipped)?;
    columns.set_item("qual", quals)?;
    columns.set_item("af", afs)?;
    columns.set_item("maf", mafs)?;
    columns.set_item("mac", macs)?;
    columns.set_item("missing_rate", missing_rates)?;
    columns.set_item("n_called", n_called)?;
    Ok(columns)
}

fn py_to_json_value(value: &Bound<'_, PyAny>) -> PyResult<serde_json::Value> {
    // PyO3 exposes bool as an int subclass, so check bool before PyInt to keep
    // filter IR values semantically stable.
    if value.is_none() {
        return Ok(serde_json::Value::Null);
    }
    if value.downcast::<PyBool>().is_ok() {
        return Ok(serde_json::Value::Bool(value.extract::<bool>()?));
    }
    if value.downcast::<PyString>().is_ok() {
        return Ok(serde_json::Value::String(value.extract::<String>()?));
    }
    if value.downcast::<PyInt>().is_ok() {
        return Ok(serde_json::Value::Number(value.extract::<i64>()?.into()));
    }
    if value.downcast::<PyFloat>().is_ok() {
        let number = serde_json::Number::from_f64(value.extract::<f64>()?).ok_or_else(|| {
            pyo3::exceptions::PyValueError::new_err("filter IR contains a non-finite float")
        })?;
        return Ok(serde_json::Value::Number(number));
    }
    if let Ok(dict) = value.downcast::<PyDict>() {
        let mut object = serde_json::Map::new();
        for (key, item) in dict.iter() {
            object.insert(key.extract::<String>()?, py_to_json_value(&item)?);
        }
        return Ok(serde_json::Value::Object(object));
    }
    if let Ok(list) = value.downcast::<PyList>() {
        return list
            .iter()
            .map(|item| py_to_json_value(&item))
            .collect::<PyResult<Vec<_>>>()
            .map(serde_json::Value::Array);
    }
    if let Ok(tuple) = value.downcast::<PyTuple>() {
        return tuple
            .iter()
            .map(|item| py_to_json_value(&item))
            .collect::<PyResult<Vec<_>>>()
            .map(serde_json::Value::Array);
    }
    Err(pyo3::exceptions::PyValueError::new_err(
        "filter IR must contain only JSON-compatible values",
    ))
}
