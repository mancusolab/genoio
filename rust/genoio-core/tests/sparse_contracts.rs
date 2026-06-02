use genoio_core::{DenseDiagnostics, SampleRecord, SparseGenotypeMatrix, VariantRecord};

fn sample(id: &str) -> SampleRecord {
    SampleRecord {
        fid: None,
        iid: id.to_string(),
        father: None,
        mother: None,
        sex: None,
        phenotype: None,
    }
}

fn variant(id: &str) -> VariantRecord {
    VariantRecord {
        chrom: "1".to_string(),
        pos: 10,
        id: id.to_string(),
        a0: "A".to_string(),
        a1: "G".to_string(),
        ref_allele: Some("A".to_string()),
        alt_allele: Some("G".to_string()),
        source_a0: "A".to_string(),
        source_a1: "G".to_string(),
        flipped: false,
        af: None,
        maf: None,
        mac: None,
        missing_rate: None,
        n_called: None,
    }
}

#[test]
fn sparse_contract_rejects_malformed_csc_arrays() {
    let samples = vec![sample("S1"), sample("S2")];
    let variants = vec![variant("rs1")];

    let bad_terminal_pointer = SparseGenotypeMatrix::new(
        2,
        1,
        vec![0, 2],
        vec![0],
        vec![1.0],
        samples.clone(),
        variants.clone(),
        DenseDiagnostics::default(),
    );
    assert!(bad_terminal_pointer
        .expect_err("terminal pointer mismatch should fail")
        .to_string()
        .contains("terminal pointer"));

    let out_of_bounds_row = SparseGenotypeMatrix::new(
        2,
        1,
        vec![0, 1],
        vec![2],
        vec![1.0],
        samples,
        variants,
        DenseDiagnostics::default(),
    );
    assert!(out_of_bounds_row
        .expect_err("row index outside n_rows should fail")
        .to_string()
        .contains("row index"));
}

#[test]
fn sparse_contract_accepts_valid_empty_columns() {
    let sparse = SparseGenotypeMatrix::new(
        2,
        2,
        vec![0, 0, 1],
        vec![1],
        vec![2.0],
        vec![sample("S1"), sample("S2")],
        vec![variant("rs1"), variant("rs2")],
        DenseDiagnostics::default(),
    )
    .expect("valid sparse matrix should pass");

    assert_eq!(sparse.n_rows, 2);
    assert_eq!(sparse.n_cols, 2);
    assert_eq!(sparse.indptr, vec![0, 0, 1]);
    assert_eq!(sparse.indices, vec![1]);
    assert_eq!(sparse.data, vec![2.0]);
}
