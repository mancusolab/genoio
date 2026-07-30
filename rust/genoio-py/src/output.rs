// pattern: Imperative Shell

use std::sync::Arc;

use arrow_array::{ArrayRef, Int64Array, RecordBatch, StringArray, StructArray};
use arrow_buffer::{BooleanBuffer, Buffer, NullBuffer, OffsetBuffer, ScalarBuffer};
use arrow_schema::{ArrowError, DataType, Field, Schema};
use genoio_core::{
    DenseLayout, MetadataOutput, NullableStringColumnBuffers, SampleMetadataBuffers,
    StringColumnBuffers, VariantMetadataBuffers,
};
use numpy::{Element, PyArray1};
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyCapsule, PyDict};
use pyo3_arrow::ffi::{to_stream_pycapsule, ArrayIterator};

use crate::errors::RustInternalError;

#[pyclass(module = "genoio._rust", skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct ArrowMetadataFrame {
    batch: RecordBatch,
}

#[pymethods]
impl ArrowMetadataFrame {
    #[pyo3(signature = (requested_schema=None))]
    fn __arrow_c_stream__<'py>(
        &self,
        py: Python<'py>,
        requested_schema: Option<Bound<'py, PyCapsule>>,
    ) -> PyResult<Bound<'py, PyCapsule>> {
        let schema = self.batch.schema();
        let fields = schema.fields().clone();
        let metadata = schema.metadata.clone();
        let arrays = vec![Ok(
            Arc::new(StructArray::from(self.batch.clone())) as ArrayRef
        )];
        let array_reader = Box::new(ArrayIterator::new(
            arrays,
            Field::new_struct("", fields, false)
                .with_metadata(metadata)
                .into(),
        ));

        to_stream_pycapsule(py, array_reader, requested_schema).map_err(PyErr::from)
    }
}

pub(crate) fn metadata_to_py(py: Python<'_>, output: MetadataOutput) -> PyResult<Py<PyDict>> {
    let MetadataOutput {
        samples,
        variants,
        capabilities,
    } = output;
    let dict = PyDict::new(py);
    dict.set_item("samples", sample_arrow_buffers_to_py(py, samples)?)?;
    dict.set_item("variants", variant_arrow_buffers_to_py(py, variants)?)?;
    dict.set_item("capabilities", source_capabilities_to_py(py, capabilities)?)?;
    Ok(dict.unbind())
}

pub(crate) fn dense_matrix_to_py(
    py: Python<'_>,
    output: genoio_core::DenseGenotypeMatrix,
    return_samples: bool,
    return_variants: bool,
) -> PyResult<Py<PyDict>> {
    let genoio_core::DenseGenotypeMatrix {
        n_samples,
        n_variants,
        values,
        layout,
        samples,
        variants,
        diagnostics: dense_diagnostics,
    } = output;

    let dict = PyDict::new(py);
    dict.set_item("values", f32_vec_to_numpy(py, values)?)?;
    dict.set_item("shape", (n_samples, n_variants))?;
    dict.set_item("values_layout", dense_layout_to_py(layout))?;
    if return_samples {
        let samples = samples.ok_or_else(|| {
            RustInternalError::new_err("dense read omitted requested sample metadata")
        })?;
        dict.set_item("samples", sample_arrow_buffers_to_py(py, samples)?)?;
    }
    if return_variants {
        let variants = variants.ok_or_else(|| {
            RustInternalError::new_err("dense read omitted requested variant metadata")
        })?;
        dict.set_item("variants", variant_arrow_buffers_to_py(py, variants)?)?;
    }

    let diagnostics = PyDict::new(py);
    diagnostics.set_item("requested_samples", dense_diagnostics.requested_samples)?;
    diagnostics.set_item("retained_samples", dense_diagnostics.retained_samples)?;
    diagnostics.set_item("missing_samples", dense_diagnostics.missing_samples)?;
    diagnostics.set_item("candidate_variants", dense_diagnostics.candidate_variants)?;
    diagnostics.set_item("retained_variants", dense_diagnostics.retained_variants)?;
    diagnostics.set_item(
        "dropped_metadata_variants",
        dense_diagnostics.dropped_metadata_variants,
    )?;
    diagnostics.set_item(
        "dropped_genotype_variants",
        dense_diagnostics.dropped_genotype_variants,
    )?;
    dict.set_item("diagnostics", diagnostics)?;

    Ok(dict.unbind())
}

pub(crate) fn sparse_matrix_to_py(
    py: Python<'_>,
    output: genoio_core::SparseGenotypeMatrix,
    return_samples: bool,
    return_variants: bool,
) -> PyResult<Py<PyDict>> {
    let dict = PyDict::new(py);
    // Core sparse buffers already use SciPy's int32 index width, so transfer
    // vector ownership directly to NumPy instead of widening through i64.
    dict.set_item("indptr", vec_to_numpy(py, output.indptr))?;
    dict.set_item("indices", vec_to_numpy(py, output.indices))?;
    dict.set_item("data", f32_vec_to_numpy(py, output.data)?)?;
    dict.set_item("shape", (output.n_rows, output.n_cols))?;
    if return_samples {
        let samples = output.samples.ok_or_else(|| {
            RustInternalError::new_err("sparse read omitted requested sample metadata")
        })?;
        dict.set_item("samples", sample_arrow_buffers_to_py(py, samples)?)?;
    }
    if return_variants {
        let variants = output.variants.ok_or_else(|| {
            RustInternalError::new_err("sparse read omitted requested variant metadata")
        })?;
        dict.set_item("variants", variant_arrow_buffers_to_py(py, variants)?)?;
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

pub(crate) fn block_output_to_py(
    py: Python<'_>,
    output: genoio_io::BlockOutput,
    return_samples: bool,
    return_variants: bool,
) -> PyResult<Py<PyDict>> {
    match output {
        genoio_io::BlockOutput::Dense(matrix) => {
            dense_matrix_to_py(py, matrix, return_samples, return_variants)
        }
        genoio_io::BlockOutput::Sparse(matrix) => {
            sparse_matrix_to_py(py, matrix, return_samples, return_variants)
        }
    }
}

fn source_capabilities_to_py(
    py: Python<'_>,
    source_capabilities: genoio_core::SourceCapabilities,
) -> PyResult<Bound<'_, PyDict>> {
    let capabilities = PyDict::new(py);
    capabilities.set_item("supports_geno", source_capabilities.supports_geno)?;
    capabilities.set_item("supports_haplo", source_capabilities.supports_haplo)?;
    capabilities.set_item("phased", source_capabilities.phased)?;
    Ok(capabilities)
}

fn dense_layout_to_py(layout: DenseLayout) -> &'static str {
    match layout {
        DenseLayout::SampleMajor => "sample_major",
        DenseLayout::VariantMajor => "variant_major",
    }
}

fn f32_vec_to_numpy(py: Python<'_>, values: Vec<f32>) -> PyResult<Bound<'_, PyAny>> {
    Ok(vec_to_numpy(py, values))
}

fn vec_to_numpy<'py, T>(py: Python<'py>, values: Vec<T>) -> Bound<'py, PyAny>
where
    T: Element,
{
    PyArray1::from_vec(py, values).into_any()
}

fn sample_arrow_buffers_to_py(
    py: Python<'_>,
    samples: SampleMetadataBuffers,
) -> PyResult<Py<ArrowMetadataFrame>> {
    Py::new(
        py,
        ArrowMetadataFrame {
            batch: sample_arrow_buffers_to_batch(samples)?,
        },
    )
}

fn variant_arrow_buffers_to_py(
    py: Python<'_>,
    variants: VariantMetadataBuffers,
) -> PyResult<Py<ArrowMetadataFrame>> {
    Py::new(
        py,
        ArrowMetadataFrame {
            batch: variant_arrow_buffers_to_batch(variants)?,
        },
    )
}

fn sample_arrow_buffers_to_batch(samples: SampleMetadataBuffers) -> PyResult<RecordBatch> {
    // Consume Rust-owned column buffers directly into Arrow arrays. This is the
    // phase-2 boundary: no SampleRecord rows are rebuilt on the Python side.
    let mut fields = vec![
        Field::new("fid", DataType::Utf8, true),
        Field::new("iid", DataType::Utf8, false),
        Field::new("father", DataType::Utf8, true),
        Field::new("mother", DataType::Utf8, true),
        Field::new("sex", DataType::Utf8, true),
        Field::new("phenotype", DataType::Utf8, true),
    ];
    let mut arrays: Vec<ArrayRef> = vec![
        Arc::new(nullable_string_array_from_buffers(samples.fids)) as ArrayRef,
        Arc::new(string_array_from_buffers(samples.iids)) as ArrayRef,
        Arc::new(nullable_string_array_from_buffers(samples.fathers)) as ArrayRef,
        Arc::new(nullable_string_array_from_buffers(samples.mothers)) as ArrayRef,
        Arc::new(nullable_string_array_from_buffers(samples.sexes)) as ArrayRef,
        Arc::new(nullable_string_array_from_buffers(samples.phenotypes)) as ArrayRef,
    ];

    if let Some(source_sample_indices) = samples.source_sample_indices {
        // Mapping columns are present only for haplotype-expanded reads, where
        // each public row needs to point back to its source sample and allele.
        fields.push(Field::new("source_sample_index", DataType::Int64, true));
        arrays.push(Arc::new(optional_usize_to_i64_array(source_sample_indices)?) as ArrayRef);
    }
    if let Some(haplotype_indices) = samples.haplotype_indices {
        fields.push(Field::new("haplotype_index", DataType::Int64, true));
        arrays.push(Arc::new(optional_usize_to_i64_array(haplotype_indices)?) as ArrayRef);
    }

    RecordBatch::try_new(Arc::new(Schema::new(fields)), arrays).map_err(arrow_error_to_py)
}

fn variant_arrow_buffers_to_batch(variants: VariantMetadataBuffers) -> PyResult<RecordBatch> {
    let schema = Arc::new(Schema::new(vec![
        Field::new("chrom", DataType::Utf8, false),
        Field::new("pos", DataType::Int64, false),
        Field::new("id", DataType::Utf8, false),
        Field::new("a0", DataType::Utf8, false),
        Field::new("a1", DataType::Utf8, false),
    ]));
    let arrays: Vec<ArrayRef> = vec![
        Arc::new(string_array_from_buffers(variants.chroms)) as ArrayRef,
        Arc::new(Int64Array::from(variants.positions)) as ArrayRef,
        Arc::new(string_array_from_buffers(variants.ids)) as ArrayRef,
        Arc::new(string_array_from_buffers(variants.a0s)) as ArrayRef,
        Arc::new(string_array_from_buffers(variants.a1s)) as ArrayRef,
    ];

    RecordBatch::try_new(schema, arrays).map_err(arrow_error_to_py)
}

fn string_array_from_buffers(column: StringColumnBuffers) -> StringArray {
    let offsets = ScalarBuffer::from(column.offsets);
    let values = Buffer::from_vec(column.values);
    // SAFETY: StringColumnBuffers is only populated via append_value(&str),
    // which appends valid UTF-8 bytes and cumulative non-negative i32 byte
    // offsets beginning at 0. There are no nulls in the public VCF variant
    // metadata columns.
    unsafe { StringArray::new_unchecked(OffsetBuffer::new_unchecked(offsets), values, None) }
}

fn nullable_string_array_from_buffers(column: NullableStringColumnBuffers) -> StringArray {
    let nulls = (!column.validity.iter().all(|&is_valid| is_valid))
        .then(|| NullBuffer::new(BooleanBuffer::from(column.validity)));
    let offsets = ScalarBuffer::from(column.offsets);
    let values = Buffer::from_vec(column.values);
    // SAFETY: NullableStringColumnBuffers uses the same append-only UTF-8 and
    // offset invariants as StringColumnBuffers, with one validity bit per row.
    unsafe { StringArray::new_unchecked(OffsetBuffer::new_unchecked(offsets), values, nulls) }
}

fn optional_usize_to_i64_array(values: Vec<Option<usize>>) -> PyResult<Int64Array> {
    let values = values
        .into_iter()
        .map(|value| {
            value
                .map(|value| {
                    i64::try_from(value).map_err(|_| {
                        pyo3::exceptions::PyOverflowError::new_err(
                            "array index exceeds supported Arrow int64 range",
                        )
                    })
                })
                .transpose()
        })
        .collect::<PyResult<Vec<_>>>()?;
    Ok(Int64Array::from(values))
}

fn arrow_error_to_py(error: ArrowError) -> PyErr {
    RustInternalError::new_err(format!("failed to build Arrow metadata payload: {error}"))
}
