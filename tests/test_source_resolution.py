import pytest


def test_vcf_missing_path_raises_invalid_source_error(tmp_path):
    import genoio

    with pytest.raises(genoio.InvalidSourceError):
        genoio.vcf(tmp_path / "missing.vcf")


def test_vcf_unsupported_extension_raises_unsupported_format_error(tmp_path):
    import genoio

    source_path = tmp_path / "cohort.txt"
    source_path.touch()

    with pytest.raises(genoio.UnsupportedFormatError):
        genoio.vcf(source_path)


def test_vcf_does_not_resolve_same_stem_plink_prefix(tmp_path):
    import genoio

    for suffix in (".bed", ".bim", ".fam"):
        (tmp_path / f"cohort{suffix}").touch()
    source_path = tmp_path / "cohort.txt"
    source_path.touch()

    with pytest.raises(genoio.UnsupportedFormatError):
        genoio.vcf(source_path)


def test_bgen_member_path_resolves_bgen_source(tmp_path):
    import genoio

    source_path = tmp_path / "cohort.bgen"
    source_path.touch()

    dataset = genoio.bgen(source_path)

    assert dataset.source.format.value == "bgen"
    assert dataset.source.path == source_path
    assert dataset.source.members == {"bgen": source_path}
    assert dataset.source.prefix == tmp_path / "cohort"


def test_bgen_prefix_resolves_optional_sample_companion(tmp_path):
    import genoio

    bgen_path = tmp_path / "cohort.bgen"
    sample_path = tmp_path / "cohort.sample"
    bgen_path.touch()
    sample_path.touch()

    dataset = genoio.bgen(tmp_path / "cohort")

    assert dataset.source.members == {"bgen": bgen_path, "sample": sample_path}


def test_bgen_dotted_member_path_preserves_full_prefix(tmp_path):
    import genoio

    source_path = tmp_path / "cohort.v1.bgen"
    source_path.touch()

    dataset = genoio.bgen(source_path)

    assert dataset.source.path == source_path
    assert dataset.source.members == {"bgen": source_path}
    assert dataset.source.prefix == tmp_path / "cohort.v1"


def test_bgen_dotted_prefix_resolves_optional_sample_companion(tmp_path):
    import genoio

    bgen_path = tmp_path / "cohort.v1.bgen"
    sample_path = tmp_path / "cohort.v1.sample"
    bgen_path.touch()
    sample_path.touch()

    dataset = genoio.bgen(tmp_path / "cohort.v1")

    assert dataset.source.path == bgen_path
    assert dataset.source.members == {"bgen": bgen_path, "sample": sample_path}
    assert dataset.source.prefix == tmp_path / "cohort.v1"


def test_bgen_missing_path_raises_invalid_source_error(tmp_path):
    import genoio

    with pytest.raises(genoio.InvalidSourceError):
        genoio.bgen(tmp_path / "missing.bgen")


def test_bgen_unsupported_extension_raises_unsupported_format_error(tmp_path):
    import genoio

    source_path = tmp_path / "cohort.vcf"
    source_path.touch()

    with pytest.raises(genoio.UnsupportedFormatError):
        genoio.bgen(source_path)


def test_bfile_missing_companion_raises_missing_companion_error(tmp_path):
    import genoio

    (tmp_path / "cohort.bed").touch()
    (tmp_path / "cohort.bim").touch()

    with pytest.raises(genoio.MissingCompanionFileError, match="cohort.fam"):
        genoio.bfile(tmp_path / "cohort")


def test_bfile_member_path_resolves_shared_prefix(tmp_path):
    import genoio

    for suffix in (".bed", ".bim", ".fam"):
        (tmp_path / f"cohort{suffix}").touch()

    dataset = genoio.bfile(tmp_path / "cohort.bed")

    assert dataset.source.format.value == "plink1"
    assert dataset.source.prefix == tmp_path / "cohort"
    assert set(dataset.source.members) == {"bed", "bim", "fam"}


def test_bfile_rejects_non_bfile_member_path(tmp_path):
    import genoio

    source_path = tmp_path / "cohort.pgen"
    source_path.touch()

    with pytest.raises(genoio.UnsupportedFormatError, match="not plink1"):
        genoio.bfile(source_path)


def test_pfile_missing_companion_raises_missing_companion_error(tmp_path):
    import genoio

    (tmp_path / "cohort.pgen").touch()
    (tmp_path / "cohort.pvar").touch()

    with pytest.raises(genoio.MissingCompanionFileError, match="cohort.psam"):
        genoio.pfile(tmp_path / "cohort")


def test_pfile_member_path_resolves_shared_prefix(tmp_path):
    import genoio

    for suffix in (".pgen", ".pvar", ".psam"):
        (tmp_path / f"cohort{suffix}").touch()

    dataset = genoio.pfile(tmp_path / "cohort.pgen")

    assert dataset.source.format.value == "plink2"
    assert dataset.source.prefix == tmp_path / "cohort"
    assert set(dataset.source.members) == {"pgen", "pvar", "psam"}


def test_pfile_resolves_compressed_pvar_when_uncompressed_is_absent(tmp_path):
    import genoio

    for suffix in (".pgen", ".pvar.zst", ".psam"):
        (tmp_path / f"cohort{suffix}").touch()

    dataset = genoio.pfile(tmp_path / "cohort")

    assert dataset.source.members["pvar"] == tmp_path / "cohort.pvar.zst"


def test_pfile_compressed_pvar_member_path_resolves_shared_prefix(tmp_path):
    import genoio

    for suffix in (".pgen", ".pvar.zst", ".psam"):
        (tmp_path / f"cohort{suffix}").touch()

    dataset = genoio.pfile(tmp_path / "cohort.pvar.zst")

    assert dataset.source.prefix == tmp_path / "cohort"
    assert dataset.source.members["pvar"] == tmp_path / "cohort.pvar.zst"


def test_pfile_prefers_uncompressed_pvar_when_both_exist(tmp_path):
    import genoio

    for suffix in (".pgen", ".pvar", ".pvar.zst", ".psam"):
        (tmp_path / f"cohort{suffix}").touch()

    dataset = genoio.pfile(tmp_path / "cohort")

    assert dataset.source.members["pvar"] == tmp_path / "cohort.pvar"


def test_pfile_rejects_non_pfile_member_path(tmp_path):
    import genoio

    source_path = tmp_path / "cohort.bed"
    source_path.touch()

    with pytest.raises(genoio.UnsupportedFormatError, match="not plink2"):
        genoio.pfile(source_path)


def test_explicit_pfile_ignores_stray_bfile_member(tmp_path):
    import genoio

    (tmp_path / "cohort.bed").touch()
    for suffix in (".pgen", ".pvar", ".psam"):
        (tmp_path / f"cohort{suffix}").touch()

    dataset = genoio.pfile(tmp_path / "cohort")

    assert dataset.source.format.value == "plink2"
    assert set(dataset.source.members) == {"pgen", "pvar", "psam"}


def test_plink2_read_rejects_invalid_empty_files(tmp_path):
    import genoio

    for suffix in (".pgen", ".pvar", ".psam"):
        (tmp_path / f"cohort{suffix}").touch()

    dataset = genoio.pfile(tmp_path / "cohort")

    with pytest.raises(genoio.InvalidSourceError):
        dataset.read()
