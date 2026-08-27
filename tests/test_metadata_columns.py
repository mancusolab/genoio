# pattern: Functional Core

import warnings

import polars as pl

from genoio._assembly import samples_frame, variants_frame


def test_samples_frame_builds_from_column_payload():
    frame = samples_frame(
        {
            "fid": ["F1", "F2"],
            "iid": ["S1", "S2"],
            "father": [None, None],
            "mother": [None, None],
            "sex": ["1", "2"],
            "phenotype": [None, None],
            "source_sample_index": [None, None],
            "haplotype_index": [None, None],
        }
    )

    assert frame.columns == ["fid", "iid", "father", "mother", "sex", "phenotype"]
    assert frame.to_dict(as_series=False) == {
        "fid": ["F1", "F2"],
        "iid": ["S1", "S2"],
        "father": [None, None],
        "mother": [None, None],
        "sex": ["1", "2"],
        "phenotype": [None, None],
    }


def test_variants_frame_builds_from_column_payload():
    frame = variants_frame(
        {
            "chrom": ["1", "1"],
            "pos": [10, 20],
            "id": ["rs1", "rs2"],
            "a0": ["A", "C"],
            "a1": ["G", "T"],
        }
    )

    assert frame.columns == ["chrom", "pos", "id", "a0", "a1"]
    assert frame.to_dict(as_series=False) == {
        "chrom": ["1", "1"],
        "pos": [10, 20],
        "id": ["rs1", "rs2"],
        "a0": ["A", "C"],
        "a1": ["G", "T"],
    }


def test_arrow_stream_metadata_does_not_emit_dimensionality_future_warning(monkeypatch):
    payload = pl.DataFrame(
        {
            "chrom": ["1"],
            "pos": [10],
            "id": ["rs1"],
            "a0": ["A"],
            "a1": ["G"],
        }
    )

    def ambiguous_arrow_conversion(stream):
        warnings.warn(
            "Arrow stream dimensionality will change in Polars 2.0",
            FutureWarning,
            stacklevel=2,
        )
        return pl.DataFrame(stream)

    monkeypatch.setattr(pl, "from_arrow", ambiguous_arrow_conversion)

    with warnings.catch_warnings():
        warnings.simplefilter("error", FutureWarning)
        frame = variants_frame(payload)

    assert frame.equals(payload)
