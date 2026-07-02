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
