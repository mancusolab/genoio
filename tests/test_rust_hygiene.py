# pattern: Imperative Shell

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

    assert "::from_matrix" not in source
    assert "fn read_dense_matrix(" not in source
    assert "fn read_sparse_matrix(" not in source


def test_rust_backends_do_not_use_migration_output_names():
    sources = "\n".join(
        path.read_text()
        for root in [
            Path("rust/genoio-core/src"),
            Path("rust/genoio-io/src"),
            Path("rust/genoio-py/src"),
        ]
        for path in root.rglob("*.rs")
    )

    assert "ArrowVariants" not in sources
    assert "arrow_variants" not in sources
    assert "with_arrow_variants" not in sources
    assert "MetadataArrowOutput" not in sources
    assert "SampleMetadataArrowBuffers" not in sources
    assert "VariantMetadataArrowBuffers" not in sources


def test_plink1_output_backend_does_not_wrap_row_matrices():
    source = Path("rust/genoio-io/src/plink/plink1.rs").read_text()

    assert "dense_matrix_to_arrow_variants" not in source
    assert "sparse_matrix_to_arrow_variants" not in source


def test_plink2_hardcall_output_backend_does_not_wrap_row_matrices():
    source = Path("rust/genoio-io/src/plink/plink2.rs").read_text()

    dense_wrapper = re.compile(
        r"pub fn read_plink2_dense_windowed[\s\S]*?"
        r"dense_matrix_to_arrow_variants"
    )
    sparse_wrapper = re.compile(
        r"pub fn read_plink2_sparse_windowed[\s\S]*?"
        r"sparse_matrix_to_arrow_variants"
    )

    assert dense_wrapper.search(source) is None
    assert sparse_wrapper.search(source) is None


def test_plink2_non_hardcall_output_backend_does_not_wrap_row_matrices():
    source = Path("rust/genoio-io/src/plink/plink2.rs").read_text()

    wrapped_entry_points = [
        (
            "read_plink2_dosage_dense_windowed",
            "dense_matrix_to_arrow_variants",
        ),
        (
            "read_plink2_haplotypes_dense_windowed",
            "dense_matrix_to_arrow_variants",
        ),
        (
            "read_plink2_haplotypes_dosage_dense_windowed",
            "dense_matrix_to_arrow_variants",
        ),
        (
            "read_plink2_haplotypes_sparse_windowed",
            "sparse_matrix_to_arrow_variants",
        ),
    ]

    for function_name, conversion_helper in wrapped_entry_points:
        wrapper = re.compile(rf"pub fn {function_name}[\s\S]*?{conversion_helper}")
        assert wrapper.search(source) is None


def test_bgen_output_backend_does_not_wrap_row_matrices():
    source = Path("rust/genoio-io/src/bgen.rs").read_text()

    wrapped_entry_points = [
        "read_bgen_dosage_dense_windowed",
        "read_bgen_haplotypes_dosage_dense_windowed",
    ]

    for function_name in wrapped_entry_points:
        wrapper = re.compile(rf"pub fn {function_name}[\s\S]*?dense_matrix_to_arrow_variants")
        assert wrapper.search(source) is None


def test_bgen_metadata_and_zstd_paths_avoid_decode_hot_spot_allocations():
    decode_source = Path("rust/genoio-io/src/bgen/decode.rs").read_text()
    header_source = Path("rust/genoio-io/src/bgen/header.rs").read_text()

    assert "zstd::stream::decode_all" not in decode_source
    assert "skip_layout2_probability_block" not in decode_source
    assert "skip_layout2_probability_block" not in header_source
    assert "skip_layout2_probability_payload_raw" in header_source


def test_bcf_output_backend_does_not_wrap_row_matrices():
    source = Path("rust/genoio-io/src/vcf.rs").read_text()

    assert "dense_matrix_to_arrow_variants" not in source
    assert "sparse_matrix_to_arrow_variants" not in source


def test_plink_metadata_paths_are_streamed():
    plink1_metadata = Path("rust/genoio-io/src/plink/plink1/metadata.rs").read_text()
    plink2_metadata = Path("rust/genoio-io/src/plink/plink2/metadata.rs").read_text()

    assert "fs::read_to_string" not in plink1_metadata
    assert "fs::read_to_string" not in plink2_metadata
    assert ".read_to_string(&mut contents)" not in plink2_metadata
    assert "data_lines = contents" not in plink2_metadata


def test_plink1_matrix_reads_do_not_materialize_full_bim_rows():
    source = Path("rust/genoio-io/src/plink/plink1.rs").read_text()

    assert "parse_bim(" not in source


def test_plink2_source_windows_do_not_synthesize_dummy_variant_records():
    dense = Path("rust/genoio-io/src/plink/plink2/dense.rs").read_text()
    sparse = Path("rust/genoio-io/src/plink/plink2/sparse.rs").read_text()

    for source in [dense, sparse]:
        assert "chrom: String::new()" not in source
        assert "id: String::new()" not in source
        assert "a0: String::new()" not in source
        assert "a1: String::new()" not in source


def test_pgen_difflist_decoder_does_not_allocate_entry_vectors():
    sources = [
        Path("rust/genoio-io/src/plink/plink2/pgen/main_track.rs").read_text(),
        Path("rust/genoio-io/src/plink/plink2/pgen/dosage_track.rs").read_text(),
    ]

    for source in sources:
        assert "let mut first_ids = Vec::" not in source
        assert "let mut entries = Vec::" not in source
        assert "Result<Vec<(usize, u8)>>" not in source
