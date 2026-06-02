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

#[pymodule]
fn _rust(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(backend_name, module)?)?;
    module.add_function(wrap_pyfunction!(read_metadata, module)?)?;
    Ok(())
}

fn member_path(members: &Bound<'_, PyDict>, key: &str) -> PyResult<PathBuf> {
    let value = members.get_item(key)?.ok_or_else(|| {
        pyo3::exceptions::PyKeyError::new_err(format!("missing source member: {key}"))
    })?;
    Ok(PathBuf::from(value.extract::<String>()?))
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
