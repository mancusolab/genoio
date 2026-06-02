import pytest


def test_missing_path_raises_invalid_source_error(tmp_path):
    import genoio

    with pytest.raises(genoio.InvalidSourceError):
        genoio.open(tmp_path / "missing.vcf")


def test_unsupported_extension_raises_unsupported_format_error(tmp_path):
    import genoio

    source_path = tmp_path / "cohort.txt"
    source_path.touch()

    with pytest.raises(genoio.UnsupportedFormatError):
        genoio.open(source_path)


def test_unsupported_extension_does_not_resolve_as_same_stem_plink_prefix(tmp_path):
    import genoio

    for suffix in (".bed", ".bim", ".fam"):
        (tmp_path / f"cohort{suffix}").touch()
    source_path = tmp_path / "cohort.txt"
    source_path.touch()

    with pytest.raises(genoio.UnsupportedFormatError):
        genoio.open(source_path)


def test_plink1_missing_companion_raises_missing_companion_error(tmp_path):
    import genoio

    (tmp_path / "cohort.bed").touch()
    (tmp_path / "cohort.bim").touch()

    with pytest.raises(genoio.MissingCompanionFileError):
        genoio.open(tmp_path / "cohort")


def test_plink1_missing_bim_raises_missing_companion_error(tmp_path):
    import genoio

    (tmp_path / "cohort.bed").touch()
    (tmp_path / "cohort.fam").touch()

    with pytest.raises(genoio.MissingCompanionFileError, match="cohort.bim"):
        genoio.open(tmp_path / "cohort")


def test_plink1_missing_fam_raises_missing_companion_error(tmp_path):
    import genoio

    (tmp_path / "cohort.bed").touch()
    (tmp_path / "cohort.bim").touch()

    with pytest.raises(genoio.MissingCompanionFileError, match="cohort.fam"):
        genoio.open(tmp_path / "cohort")


def test_plink2_missing_companion_raises_missing_companion_error(tmp_path):
    import genoio

    (tmp_path / "cohort.pgen").touch()
    (tmp_path / "cohort.pvar").touch()

    with pytest.raises(genoio.MissingCompanionFileError):
        genoio.open(tmp_path / "cohort")


def test_prefix_matching_plink1_and_plink2_raises_ambiguous_source_error(tmp_path):
    import genoio

    for suffix in (".bed", ".bim", ".fam", ".pgen", ".pvar", ".psam"):
        (tmp_path / f"cohort{suffix}").touch()

    with pytest.raises(genoio.AmbiguousSourceError):
        genoio.open(tmp_path / "cohort")


def test_resolves_complete_plink1_prefix(tmp_path):
    import genoio

    for suffix in (".bed", ".bim", ".fam"):
        (tmp_path / f"cohort{suffix}").touch()

    dataset = genoio.open(tmp_path / "cohort")

    assert dataset.source.format.value == "plink1"
    assert set(dataset.source.members) == {"bed", "bim", "fam"}


def test_resolves_complete_plink2_member_path(tmp_path):
    import genoio

    for suffix in (".pgen", ".pvar", ".psam"):
        (tmp_path / f"cohort{suffix}").touch()

    dataset = genoio.open(tmp_path / "cohort.pgen")

    assert dataset.source.format.value == "plink2"
    assert dataset.source.prefix == tmp_path / "cohort"


def test_resolves_complete_plink2_prefix(tmp_path):
    import genoio

    for suffix in (".pgen", ".pvar", ".psam"):
        (tmp_path / f"cohort{suffix}").touch()

    dataset = genoio.open(tmp_path / "cohort")

    assert dataset.source.format.value == "plink2"
    assert set(dataset.source.members) == {"pgen", "pvar", "psam"}


def test_plink2_read_raises_deferred_decode_error(tmp_path):
    import genoio

    for suffix in (".pgen", ".pvar", ".psam"):
        (tmp_path / f"cohort{suffix}").touch()

    dataset = genoio.open(tmp_path / "cohort")

    with pytest.raises(genoio.UnsupportedFormatError, match="PLINK2 decode is deferred"):
        dataset.read()


def test_explicit_plink1_format_on_non_plink_source_raises_unsupported_format_error(tmp_path):
    import genoio

    source_path = tmp_path / "cohort.vcf"
    source_path.touch()

    with pytest.raises(genoio.UnsupportedFormatError):
        genoio.open(source_path, format="plink1")
