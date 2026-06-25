use genoio_core::{
    append_sparse_column, DenseDiagnostics, SampleRecord, SparseGenotypeMatrix,
    SparseGenotypeMatrixArrowVariants, VariantRecord,
};

fn sample(id: &str) -> SampleRecord {
    SampleRecord {
        fid: None,
        iid: id.to_string(),
        father: None,
        mother: None,
        sex: None,
        phenotype: None,
        source_sample_index: None,
        haplotype_index: None,
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
        qual: None,
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
    let indptr: Vec<i32> = vec![0, 0, 1];
    let indices: Vec<i32> = vec![1];
    let sparse = SparseGenotypeMatrix::new(
        2,
        2,
        indptr,
        indices,
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

#[test]
fn sparse_arrow_contract_rejects_negative_indices() {
    let negative_pointer = SparseGenotypeMatrixArrowVariants::new(
        2,
        1,
        vec![0, -1],
        vec![],
        vec![],
        None,
        None,
        DenseDiagnostics::default(),
    );
    assert!(negative_pointer
        .expect_err("negative indptr should fail")
        .to_string()
        .contains("must be nonnegative"));

    let negative_row = SparseGenotypeMatrixArrowVariants::new(
        2,
        1,
        vec![0, 1],
        vec![-1],
        vec![1.0],
        None,
        None,
        DenseDiagnostics::default(),
    );
    assert!(negative_row
        .expect_err("negative row index should fail")
        .to_string()
        .contains("must be nonnegative"));
}

#[test]
fn sparse_arrow_contract_rejects_dimensions_outside_i32_index_range() {
    let too_many_rows = SparseGenotypeMatrixArrowVariants::new(
        i32::MAX as usize + 1,
        0,
        vec![0],
        vec![],
        vec![],
        None,
        None,
        DenseDiagnostics::default(),
    );
    assert!(too_many_rows
        .expect_err("n_rows outside i32 range should fail")
        .to_string()
        .contains("exceeds sparse int32 index range"));

    let too_many_columns = SparseGenotypeMatrixArrowVariants::new(
        0,
        i32::MAX as usize,
        vec![0],
        vec![],
        vec![],
        None,
        None,
        DenseDiagnostics::default(),
    );
    assert!(too_many_columns
        .expect_err("n_cols + 1 outside i32 range should fail")
        .to_string()
        .contains("exceeds sparse int32 index range"));
}

#[test]
fn append_sparse_column_emits_i32_indices() {
    let mut indptr: Vec<i32> = vec![0];
    let mut indices: Vec<i32> = Vec::new();
    let mut data = Vec::new();

    append_sparse_column(&mut indptr, &mut indices, &mut data, &[0.0, 1.0, 2.0])
        .expect("small sparse column should append");

    assert_eq!(indptr, vec![0, 2]);
    assert_eq!(indices, vec![1, 2]);
    assert_eq!(data, vec![1.0, 2.0]);
}
