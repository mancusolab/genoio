import re
from pathlib import Path


def test_internal_contract_errors_map_to_private_internal_error():
    source = Path("rust/genoio-py/src/lib.rs").read_text()
    branch_match = re.search(r"GenoioError::InternalContract\s*\{ \.\. \}\s*=>.*", source)

    assert branch_match is not None
    branch = branch_match.group(0)
    assert "RustInternalError::new_err" in branch
    assert "PyRuntimeError" not in branch


def test_text_vcf_backend_has_no_row_variant_sink():
    text_sources = [
        Path("rust/genoio-io/src/vcf/text.rs").read_text(),
        Path("rust/genoio-io/src/vcf/text/sparse.rs").read_text(),
    ]
    combined = "\n".join(text_sources)

    assert "VariantMetadataSinkKind::Records" not in combined
    assert "TextDenseReadOutput::Records" not in combined
    assert "TextSparseReadOutput::Records" not in combined
    assert 'into_records("sparse")' not in combined


def test_pyo3_adapter_does_not_normalize_row_matrices():
    source = Path("rust/genoio-py/src/lib.rs").read_text()

    assert "DenseGenotypeMatrixArrowVariants::from_matrix" not in source
    assert "SparseGenotypeMatrixArrowVariants::from_matrix" not in source
    assert "fn read_dense_matrix(" not in source
    assert "fn read_sparse_matrix(" not in source


def test_plink1_arrow_backend_does_not_wrap_row_matrices():
    source = Path("rust/genoio-io/src/plink/plink1.rs").read_text()

    assert "dense_matrix_to_arrow_variants" not in source
    assert "sparse_matrix_to_arrow_variants" not in source
