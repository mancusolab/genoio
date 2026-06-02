// pattern: Imperative Shell

use std::path::PathBuf;

use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};

#[pyfunction]
fn backend_name() -> &'static str {
    genoio_io::backend_name()
}

#[pyfunction]
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
fn read_dense(
    py: Python<'_>,
    format: &str,
    members: &Bound<'_, PyDict>,
    options: &Bound<'_, PyDict>,
) -> PyResult<Py<PyDict>> {
    let requested_samples = samples_option(options)?;
    let output = match format {
        "vcf" | "bcf" => {
            let key = if format == "vcf" { "vcf" } else { "bcf" };
            let path = member_path(members, key)?;
            genoio_io::read_vcf_dense(&path, requested_samples.as_deref())
        }
        "plink1" => {
            let bed = member_path(members, "bed")?;
            let bim = member_path(members, "bim")?;
            let fam = member_path(members, "fam")?;
            genoio_io::read_plink1_dense(&bed, &bim, &fam, requested_samples.as_deref())
        }
        other => {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "unsupported dense format: {other}"
            )));
        }
    }
    .map_err(|error| pyo3::exceptions::PyValueError::new_err(error.to_string()))?;

    dense_to_py(py, output)
}

#[pymodule]
fn _rust(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(backend_name, module)?)?;
    module.add_function(wrap_pyfunction!(read_metadata, module)?)?;
    module.add_function(wrap_pyfunction!(read_dense, module)?)?;
    Ok(())
}

fn member_path(members: &Bound<'_, PyDict>, key: &str) -> PyResult<PathBuf> {
    let value = members.get_item(key)?.ok_or_else(|| {
        pyo3::exceptions::PyKeyError::new_err(format!("missing source member: {key}"))
    })?;
    Ok(PathBuf::from(value.extract::<String>()?))
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

fn metadata_to_py(py: Python<'_>, output: genoio_core::MetadataOutput) -> PyResult<Py<PyDict>> {
    let dict = PyDict::new(py);
    let samples = PyList::empty(py);
    for sample in output.samples {
        let sample_dict = PyDict::new(py);
        sample_dict.set_item("fid", sample.fid)?;
        sample_dict.set_item("iid", sample.iid)?;
        sample_dict.set_item("father", sample.father)?;
        sample_dict.set_item("mother", sample.mother)?;
        sample_dict.set_item("sex", sample.sex)?;
        sample_dict.set_item("phenotype", sample.phenotype)?;
        samples.append(sample_dict)?;
    }

    let variants = PyList::empty(py);
    for variant in output.variants {
        let variant_dict = PyDict::new(py);
        variant_dict.set_item("chrom", variant.chrom)?;
        variant_dict.set_item("pos", variant.pos)?;
        variant_dict.set_item("id", variant.id)?;
        variant_dict.set_item("a0", variant.a0)?;
        variant_dict.set_item("a1", variant.a1)?;
        variant_dict.set_item("ref_allele", variant.ref_allele)?;
        variant_dict.set_item("alt_allele", variant.alt_allele)?;
        variant_dict.set_item("source_a0", variant.source_a0)?;
        variant_dict.set_item("source_a1", variant.source_a1)?;
        variant_dict.set_item("flipped", variant.flipped)?;
        variants.append(variant_dict)?;
    }

    let capabilities = PyDict::new(py);
    capabilities.set_item("supports_geno", output.capabilities.supports_geno)?;
    capabilities.set_item("supports_haplo", output.capabilities.supports_haplo)?;
    capabilities.set_item("phased", output.capabilities.phased)?;

    dict.set_item("samples", samples)?;
    dict.set_item("variants", variants)?;
    dict.set_item("capabilities", capabilities)?;
    Ok(dict.unbind())
}

fn dense_to_py(py: Python<'_>, output: genoio_core::DenseGenotypeMatrix) -> PyResult<Py<PyDict>> {
    let dict = PyDict::new(py);
    dict.set_item("values", output.values)?;
    dict.set_item("shape", (output.n_samples, output.n_variants))?;
    dict.set_item("missing_mask", output.missing_mask)?;
    dict.set_item("samples", sample_records_to_py(py, output.samples)?)?;
    dict.set_item("variants", variant_records_to_py(py, output.variants)?)?;

    let diagnostics = PyDict::new(py);
    diagnostics.set_item("requested_samples", output.diagnostics.requested_samples)?;
    diagnostics.set_item("retained_samples", output.diagnostics.retained_samples)?;
    diagnostics.set_item("missing_samples", output.diagnostics.missing_samples)?;
    dict.set_item("diagnostics", diagnostics)?;

    Ok(dict.unbind())
}

fn sample_records_to_py(
    py: Python<'_>,
    samples: Vec<genoio_core::SampleRecord>,
) -> PyResult<Bound<'_, PyList>> {
    let py_samples = PyList::empty(py);
    for sample in samples {
        let sample_dict = PyDict::new(py);
        sample_dict.set_item("fid", sample.fid)?;
        sample_dict.set_item("iid", sample.iid)?;
        sample_dict.set_item("father", sample.father)?;
        sample_dict.set_item("mother", sample.mother)?;
        sample_dict.set_item("sex", sample.sex)?;
        sample_dict.set_item("phenotype", sample.phenotype)?;
        py_samples.append(sample_dict)?;
    }
    Ok(py_samples)
}

fn variant_records_to_py(
    py: Python<'_>,
    variants: Vec<genoio_core::VariantRecord>,
) -> PyResult<Bound<'_, PyList>> {
    let py_variants = PyList::empty(py);
    for variant in variants {
        let variant_dict = PyDict::new(py);
        variant_dict.set_item("chrom", variant.chrom)?;
        variant_dict.set_item("pos", variant.pos)?;
        variant_dict.set_item("id", variant.id)?;
        variant_dict.set_item("a0", variant.a0)?;
        variant_dict.set_item("a1", variant.a1)?;
        variant_dict.set_item("ref_allele", variant.ref_allele)?;
        variant_dict.set_item("alt_allele", variant.alt_allele)?;
        variant_dict.set_item("source_a0", variant.source_a0)?;
        variant_dict.set_item("source_a1", variant.source_a1)?;
        variant_dict.set_item("flipped", variant.flipped)?;
        py_variants.append(variant_dict)?;
    }
    Ok(py_variants)
}
