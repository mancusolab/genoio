# pattern: Imperative Shell

import ast
import sys
from importlib.metadata import metadata, version
from importlib.resources import files
from pathlib import Path

import numpy as np
import pytest
from fixture_writers import (
    _bgen_sample_identifier_block,
    _bgen_variant_identifying_data,
    write_bgen_dosage,
    write_fixed_width_phased_dosage_plink2,
    write_fixed_width_plink2,
    write_phased_dosage_plink2,
    write_phased_hardcall_plink2,
)

FIXTURE_ROOT = Path(__file__).parent / "fixtures"


def write_blocks_vcf(tmp_path: Path) -> Path:
    path = tmp_path / "blocks.vcf"
    path.write_text(
        """\
##fileformat=VCFv4.2
##contig=<ID=1>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
1\t10\trs1\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1
1\t20\trs2\tC\tT\t.\tPASS\t.\tGT\t0/1\t1/1
"""
    )
    return path


def placeholder_bgen_dataset(tmp_path: Path):
    import genoio

    path = tmp_path / "cohort.bgen"
    path.touch()
    return genoio.bgen(path)


def write_invalid_phase_probability_bgen(tmp_path: Path) -> Path:
    path = tmp_path / "invalid_phase.bgen"
    contents = bytearray()
    flags = (2 << 2) | (1 << 31)
    contents.extend((20).to_bytes(4, "little"))
    contents.extend((20).to_bytes(4, "little"))
    contents.extend((1).to_bytes(4, "little"))
    contents.extend((2).to_bytes(4, "little"))
    contents.extend(b"bgen")
    contents.extend(flags.to_bytes(4, "little"))
    contents.extend(_bgen_sample_identifier_block(["sample_1", "sample_2"]))
    variant_offset = len(contents) - 4
    contents[0:4] = variant_offset.to_bytes(4, "little")
    contents.extend(_bgen_variant_identifying_data("var1", "rs1", "1", 10, ["A", "G"]))
    probability_payload = bytearray()
    probability_payload.extend((2).to_bytes(4, "little"))
    probability_payload.extend((2).to_bytes(2, "little"))
    probability_payload.extend((2).to_bytes(1, "little"))
    probability_payload.extend((2).to_bytes(1, "little"))
    probability_payload.extend([2, 2])
    probability_payload.extend((2).to_bytes(1, "little"))
    probability_payload.extend((8).to_bytes(1, "little"))
    contents.extend(len(probability_payload).to_bytes(4, "little"))
    contents.extend(probability_payload)
    path.write_bytes(contents)
    return path


def write_layout1_bgen(tmp_path: Path) -> Path:
    path = tmp_path / "layout1.bgen"
    contents = bytearray()
    flags = (1 << 2) | (1 << 31)
    contents.extend((20).to_bytes(4, "little"))
    contents.extend((20).to_bytes(4, "little"))
    contents.extend((0).to_bytes(4, "little"))
    contents.extend((1).to_bytes(4, "little"))
    contents.extend(b"bgen")
    contents.extend(flags.to_bytes(4, "little"))
    contents.extend(_bgen_sample_identifier_block(["sample_1"]))
    variant_offset = len(contents) - 4
    contents[0:4] = variant_offset.to_bytes(4, "little")
    path.write_bytes(contents)
    return path


def test_import_exposes_public_names_without_reference_packages():
    import genoio

    expected_names = {
        "__version__",
        "vcf",
        "bfile",
        "bgen",
        "pfile",
        "chrom",
        "region",
        "snp",
        "biallelic",
        "maf",
        "qual",
        "mac",
        "missing_rate",
        "polymorphic",
        "id_in",
        "GenoioError",
        "SourceResolutionError",
        "MissingCompanionFileError",
        "UnsupportedFormatError",
        "InvalidSourceError",
        "UnsupportedRepresentation",
        "InvalidOptionError",
        "InternalError",
    }

    assert expected_names <= set(dir(genoio))
    assert "jaxqtl" not in sys.modules
    assert "linear_dag" not in sys.modules


def test_python_version_export_matches_installed_metadata():
    import genoio

    assert genoio.__version__ == version("genoio")


def test_package_declares_inline_type_information():
    assert files("genoio").joinpath("py.typed").is_file()


def test_native_stub_matches_runtime_exports():
    from genoio import _rust

    stub_module = ast.parse(Path("src/genoio/_rust.pyi").read_text())
    stub_exports = {
        node.name
        for node in stub_module.body
        if isinstance(node, ast.ClassDef | ast.FunctionDef) and not node.name.startswith("_")
    }
    runtime_exports = {name for name in dir(_rust) if not name.startswith("_")}

    assert stub_exports == runtime_exports


def test_package_metadata_includes_release_urls_and_license():
    package_metadata = metadata("genoio")
    project_urls = package_metadata.get_all("Project-URL") or []

    assert package_metadata["License-File"] == "LICENSE"
    assert "Documentation, https://mancusolab.github.io/genoio" in project_urls
    assert "Repository, https://github.com/mancusolab/genoio" in project_urls


def test_open_returns_lightweight_dataset_for_vcf(tmp_path):
    import genoio

    source_path = tmp_path / "cohort.vcf"
    source_path.touch()

    dataset = genoio.vcf(source_path)

    assert dataset.source.path == source_path
    assert dataset.source.format.value == "vcf"


def test_dataset_read_recognizes_sparse_options_and_validates_missing_policy(tmp_path):
    import genoio

    source_path = tmp_path / "cohort.vcf.gz"
    source_path.touch()
    dataset = genoio.vcf(source_path)

    with pytest.raises(genoio.InvalidOptionError, match="sparse missing values"):
        dataset.read(sparse="csc", missing="nan")


def test_dataset_read_accepts_explicit_hardcall_dosage_source(tmp_path):
    import genoio

    dataset = genoio.vcf(write_blocks_vcf(tmp_path))

    np.testing.assert_array_equal(dataset.read(), dataset.read(dosage="hardcall"))


def test_dataset_read_rejects_unknown_dosage_source(tmp_path):
    import genoio

    dataset = genoio.vcf(write_blocks_vcf(tmp_path))

    with pytest.raises(genoio.InvalidOptionError, match="unsupported dosage source"):
        dataset.read(dosage="posterior")


def test_dataset_read_rejects_dosage_source_for_formats_without_reader_support(tmp_path):
    import genoio

    dataset = genoio.bfile(FIXTURE_ROOT / "plink1" / "tiny")

    with pytest.raises(genoio.UnsupportedRepresentation, match="does not support dosage-backed genotype reads"):
        dataset.read(dosage="dosage")


def test_vcf_dataset_read_rejects_haplotype_dosage_source(tmp_path):
    import genoio

    dataset = genoio.vcf(write_blocks_vcf(tmp_path))

    with pytest.raises(genoio.UnsupportedRepresentation, match="VCF haplotype dosage.*hardcall GT"):
        dataset.read(kind="haplo", dosage="dosage")


@pytest.mark.parametrize(
    ("dataset_factory", "read_options", "expected_shape", "expected_first_column"),
    [
        pytest.param(
            lambda genoio, tmp_path: genoio.pfile(write_phased_hardcall_plink2(tmp_path)),
            {"kind": "haplo"},
            (6, 2),
            None,
            id="plink2-hardcall",
        ),
        pytest.param(
            lambda genoio, tmp_path: genoio.pfile(write_phased_dosage_plink2(tmp_path)),
            {"kind": "haplo", "dosage": "dosage"},
            (6, 2),
            None,
            id="plink2-dosage",
        ),
        pytest.param(
            lambda genoio, tmp_path: genoio.pfile(write_fixed_width_phased_dosage_plink2(tmp_path)),
            {"kind": "haplo", "dosage": "dosage"},
            (6, 2),
            [0.25, 0.75, 0.0, 0.5, 1.0, 1.0],
            id="plink2-fixed-width-phased-dosage",
        ),
        pytest.param(
            lambda genoio, tmp_path: genoio.bgen(write_bgen_dosage(tmp_path, phased=True)),
            {"kind": "haplo", "dosage": "dosage"},
            (4, 2),
            None,
            id="bgen-phased-dosage",
        ),
    ],
)
def test_dataset_read_haplotype_backends_reach_supported_paths(
    tmp_path,
    dataset_factory,
    read_options,
    expected_shape,
    expected_first_column,
):
    import genoio

    dataset = dataset_factory(genoio, tmp_path)

    H = dataset.read(**read_options)

    assert H.shape == expected_shape
    if expected_first_column is not None:
        np.testing.assert_allclose(H[:, 0], expected_first_column, atol=2.0 / 32768.0)


@pytest.mark.parametrize(
    ("read_options", "match"),
    [
        pytest.param({"dosage": "hardcall"}, "hardcall", id="hardcall-genotype"),
        pytest.param({"sparse": True, "dosage": "hardcall"}, "sparse", id="sparse-hardcall-genotype"),
        pytest.param({"kind": "haplo"}, "hardcall haplotype", id="hardcall-haplotype"),
    ],
)
def test_bgen_dataset_read_rejects_unsupported_hardcall_modes(tmp_path, read_options, match):
    import genoio

    dataset = placeholder_bgen_dataset(tmp_path)

    with pytest.raises(genoio.UnsupportedRepresentation, match=match):
        dataset.read(**read_options)


def test_bgen_dataset_read_default_haplotype_does_not_imply_dosage(tmp_path):
    import genoio

    dataset = genoio.bgen(write_bgen_dosage(tmp_path, phased=True))

    with pytest.raises(genoio.UnsupportedRepresentation, match="hardcall haplotype"):
        dataset.read(kind="haplo")


@pytest.mark.parametrize(
    "read_options",
    [
        pytest.param({"kind": "haplo"}, id="dense"),
        pytest.param({"kind": "haplo", "sparse": True}, id="sparse"),
    ],
)
def test_plink2_dataset_read_rejects_fixed_width_hardcall_haplotypes(tmp_path, read_options):
    import genoio

    dataset = genoio.pfile(write_fixed_width_plink2(tmp_path))

    with pytest.raises(genoio.UnsupportedRepresentation, match="variable-width explicit phased records"):
        dataset.read(**read_options)


def test_plink2_dataset_read_phased_hardcall_as_haplotype_dosage_is_unsupported(tmp_path):
    import genoio

    dataset = genoio.pfile(write_phased_hardcall_plink2(tmp_path))

    with pytest.raises(genoio.UnsupportedRepresentation, match="explicit phased dosage values"):
        dataset.read(kind="haplo", dosage="dosage")


def test_bgen_dataset_read_rejects_sparse_haplotype_dosage_with_dense_mode_message(tmp_path):
    import genoio

    dataset = genoio.bgen(write_bgen_dosage(tmp_path, phased=True))

    with pytest.raises(genoio.UnsupportedRepresentation, match="sparse haplotype reads.*dense"):
        dataset.read(kind="haplo", dosage="dosage", sparse=True)


def test_bgen_dataset_read_dosage_rejects_invalid_placeholder_source(tmp_path):
    import genoio

    dataset = placeholder_bgen_dataset(tmp_path)

    with pytest.raises(genoio.InvalidSourceError, match="bgen"):
        dataset.read(dosage="dosage")


def test_bgen_read_maps_unsupported_probability_representation(tmp_path):
    import genoio

    dataset = genoio.bgen(write_invalid_phase_probability_bgen(tmp_path))

    with pytest.raises(genoio.UnsupportedRepresentation, match="phased probability value"):
        dataset.read(dosage="dosage")


def test_bgen_metadata_skips_unsupported_probability_representation(tmp_path):
    import genoio

    dataset = genoio.bgen(write_invalid_phase_probability_bgen(tmp_path))

    variants = dataset.variants()

    assert variants["id"].to_list() == ["rs1"]


def test_bgen_dataset_metadata_maps_unsupported_layout(tmp_path):
    import genoio

    dataset = genoio.bgen(write_layout1_bgen(tmp_path))

    with pytest.raises(genoio.UnsupportedRepresentation, match="layout"):
        dataset.variants()


def test_dataset_blocks_accepts_explicit_hardcall_dosage_source(tmp_path):
    import genoio

    dataset = genoio.vcf(write_blocks_vcf(tmp_path))

    blocks = list(dataset.iter_blocks(1, dosage="hardcall"))

    assert len(blocks) == 2
    np.testing.assert_array_equal(blocks[0], np.array([[0.0], [1.0]], dtype=np.float32))


def test_matrix_only_read_does_not_assemble_metadata_frames(monkeypatch, tmp_path):
    import genoio
    import genoio._api as api

    source_path = tmp_path / "cohort.vcf"
    source_path.touch()
    dataset = genoio.vcf(source_path)

    def fake_read_dense(format, members, options):
        assert options["return_samples"] is False
        assert options["return_variants"] is False
        assert options["missing"] == "nan"
        return {
            "values": [0.0, 1.0],
            "shape": (1, 2),
            "diagnostics": {},
        }

    def fail_frame_assembly(records):
        raise AssertionError("metadata frames should not be assembled for matrix-only reads")

    monkeypatch.setattr(api._rust, "read_dense", fake_read_dense)
    monkeypatch.setattr(api, "samples_frame", fail_frame_assembly)
    monkeypatch.setattr(api, "variants_frame", fail_frame_assembly)

    observed = dataset.read()

    np.testing.assert_array_equal(observed, np.array([[0.0, 1.0]], dtype=np.float32))


def test_dataset_read_validates_source_support_with_rust(monkeypatch, tmp_path):
    import genoio
    import genoio._api as api

    calls = []

    def fake_validate_read_support(format, kind, dosage, sparse):
        calls.append((format, kind, dosage, sparse))

    def fake_read_dense(format, members, options):
        assert options["missing"] == "nan"
        return {
            "values": np.array([0.0], dtype=np.float32),
            "shape": (1, 1),
        }

    monkeypatch.setattr(api._rust, "validate_read_support", fake_validate_read_support)
    monkeypatch.setattr(api._rust, "read_dense", fake_read_dense)
    dataset = genoio.vcf(write_blocks_vcf(tmp_path))

    dataset.read()

    assert calls == [("vcf", "geno", "hardcall", False)]


def test_dataset_read_maps_rust_support_validation_errors(monkeypatch, tmp_path):
    import genoio
    import genoio._api as api
    from genoio import _rust

    def fake_validate_read_support(format, kind, dosage, sparse):
        raise _rust.RustUnsupportedRepresentationError("unsupported by Rust policy")

    monkeypatch.setattr(api._rust, "validate_read_support", fake_validate_read_support)
    dataset = genoio.vcf(write_blocks_vcf(tmp_path))

    with pytest.raises(genoio.UnsupportedRepresentation, match="unsupported by Rust policy"):
        dataset.read()


def test_sparse_dosage_support_policy_comes_from_rust(monkeypatch, tmp_path):
    import genoio
    import genoio._api as api
    from genoio import _rust

    def fake_validate_read_support(format, kind, dosage, sparse):
        raise _rust.RustUnsupportedRepresentationError("rust sparse dosage policy")

    monkeypatch.setattr(api._rust, "validate_read_support", fake_validate_read_support)
    dataset = genoio.vcf(write_blocks_vcf(tmp_path))

    with pytest.raises(genoio.UnsupportedRepresentation, match="rust sparse dosage policy"):
        dataset.read(dosage="dosage", sparse=True)


def test_rust_dense_read_returns_numpy_buffers(tmp_path):
    import genoio
    import genoio._api as api

    path = write_blocks_vcf(tmp_path)
    dataset = genoio.vcf(path)
    members = {key: str(path) for key, path in dataset.source.members.items()}

    result = api._rust.read_dense(
        dataset.source.format.value,
        members,
        {
            "samples": None,
            "variants": None,
            "variant_window": None,
            "dosage": "hardcall",
            "missing": "raise",
            "return_samples": False,
            "return_variants": False,
            "matrix_only": True,
        },
    )

    assert isinstance(result["values"], np.ndarray)
    assert result["values"].dtype == np.float32
    assert "missing_mask" not in result
    assert "missing_indices" not in result

    metadata_result = api._rust.read_dense(
        dataset.source.format.value,
        members,
        {
            "samples": None,
            "variants": None,
            "variant_window": None,
            "dosage": "hardcall",
            "missing": "raise",
            "return_samples": True,
            "return_variants": True,
            "matrix_only": True,
        },
    )

    assert metadata_result["shape"] == (2, 2)
    assert api.samples_frame(metadata_result["samples"])["iid"].to_list() == ["S1", "S2"]
    assert api.variants_frame(metadata_result["variants"])["id"].to_list() == ["rs1", "rs2"]


def test_dataset_blocks_accepts_read_options_and_validates_size(tmp_path):
    import genoio

    source_path = tmp_path / "cohort.vcf.gz"
    source_path.touch()
    dataset = genoio.vcf(source_path)

    read_options = {
        "kind": "geno",
        "sparse": False,
        "variants": None,
        "samples": ("s2", "s1"),
        "missing": "raise",
        "dtype": "uint8",
        "return_samples": False,
        "return_variants": True,
    }

    block_iterator = dataset.iter_blocks(8192, **read_options)

    assert iter(block_iterator) is block_iterator

    with pytest.raises(genoio.InvalidOptionError, match="block size"):
        dataset.iter_blocks(0, **read_options)


def test_dataset_variants_accepts_documented_default_stats_keyword():
    import genoio

    dataset = genoio.vcf(FIXTURE_ROOT / "vcf" / "tiny.vcf")

    variants = dataset.variants(stats=None)

    assert variants["id"].to_list() == ["rs1", "rs2", "indel1"]


def test_dataset_variants_rejects_stats_until_stat_metadata_is_implemented(tmp_path):
    import genoio

    source_path = tmp_path / "cohort.vcf.gz"
    source_path.touch()
    dataset = genoio.vcf(source_path)

    with pytest.raises(genoio.InvalidOptionError, match="variant stats"):
        dataset.variants(stats=["maf"])


@pytest.mark.parametrize(
    ("read_options", "error_type"),
    [
        pytest.param({"kind": "unsupported"}, "UnsupportedRepresentation", id="unsupported-kind"),
        pytest.param({"sparse": "unsupported"}, "InvalidOptionError", id="unsupported-sparse"),
    ],
)
def test_dataset_read_rejects_unsupported_representation_options(tmp_path, read_options, error_type):
    import genoio

    source_path = tmp_path / "cohort.vcf.gz"
    source_path.touch()
    dataset = genoio.vcf(source_path)

    with pytest.raises(getattr(genoio, error_type)):
        dataset.read(**read_options)


def test_private_rust_internal_error_maps_to_public_internal_error():
    import genoio
    from genoio import _api, _rust

    public_error = _api._public_rust_error(_rust.RustInternalError("panic detail"))

    assert isinstance(public_error, genoio.InternalError)
    assert str(public_error) == "panic detail"


def test_filter_helpers_build_serializable_expressions():
    import genoio

    expression = genoio.chrom("1") & genoio.region("1:10-20") & ~genoio.missing_rate(max=0.05)

    assert expression.to_ir() == {
        "op": "and",
        "left": {
            "op": "and",
            "left": {"op": "predicate", "name": "chrom", "params": {"value": "1"}},
            "right": {"op": "predicate", "name": "region", "params": {"value": "1:10-20"}},
        },
        "right": {
            "op": "not",
            "expr": {"op": "predicate", "name": "missing_rate", "params": {"max": 0.05}},
        },
    }


def test_region_rejects_malformed_region_syntax():
    import genoio

    with pytest.raises(genoio.InvalidOptionError):
        genoio.region("not-a-region")


def test_snp_helper_is_zero_argument_snp_only_predicate():
    import genoio

    assert genoio.snp().to_ir() == {"op": "predicate", "name": "snp", "params": {}}


def test_id_in_helper_is_variant_id_matching_predicate():
    import genoio

    assert genoio.id_in(["rs2", "rs1"]).to_ir() == {
        "op": "predicate",
        "name": "id_in",
        "params": {"values": ["rs2", "rs1"]},
    }
