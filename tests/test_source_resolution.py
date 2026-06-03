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
