fn variant(id: &str, chrom: &str, pos: u32, a0: &str, a1: &str) -> genoio_core::VariantRecord {
    genoio_core::VariantRecord {
        chrom: chrom.to_string(),
        pos,
        id: id.to_string(),
        a0: a0.to_string(),
        a1: a1.to_string(),
        ref_allele: Some(a0.to_string()),
        alt_allele: Some(a1.to_string()),
        source_a0: a0.to_string(),
        source_a1: a1.to_string(),
        flipped: false,
        af: None,
        maf: None,
        mac: None,
        missing_rate: None,
        n_called: None,
    }
}

#[test]
fn filter_ir_deserializes_composed_predicates() {
    let filter = genoio_core::VariantFilter::from_json_value(serde_json::json!({
        "op": "and",
        "left": {"op": "predicate", "name": "chrom", "params": {"value": "1"}},
        "right": {
            "op": "not",
            "expr": {"op": "predicate", "name": "id_in", "params": {"values": ["rs2"]}}
        }
    }))
    .expect("filter IR should deserialize");

    assert_eq!(
        filter.metadata_decision(&variant("rs1", "1", 10, "A", "G")),
        Some(true)
    );
    assert_eq!(
        filter.metadata_decision(&variant("rs2", "1", 20, "C", "T")),
        Some(false)
    );
    assert_eq!(
        filter.metadata_decision(&variant("rs3", "2", 30, "A", "G")),
        Some(false)
    );
}

#[test]
fn filter_ir_rejects_unknown_predicate_names_and_invalid_shapes() {
    let unknown = genoio_core::VariantFilter::from_json_value(serde_json::json!({
        "op": "predicate",
        "name": "python_callback",
        "params": {}
    }));
    assert!(unknown
        .expect_err("unknown predicate should fail")
        .to_string()
        .contains("unknown predicate"));

    let invalid_region = genoio_core::VariantFilter::from_json_value(serde_json::json!({
        "op": "predicate",
        "name": "region",
        "params": {"value": "1:20-10"}
    }));
    assert!(invalid_region
        .expect_err("invalid region should fail")
        .to_string()
        .contains("region"));
}

#[test]
fn genotype_predicates_are_not_used_as_metadata_drop_decisions() {
    let filter = genoio_core::VariantFilter::from_json_value(serde_json::json!({
        "op": "and",
        "left": {"op": "predicate", "name": "chrom", "params": {"value": "1"}},
        "right": {"op": "predicate", "name": "maf", "params": {"min": 0.2}}
    }))
    .expect("filter IR should deserialize");

    assert!(filter.requires_genotype_stats());
    assert_eq!(
        filter.metadata_decision(&variant("rs1", "1", 10, "A", "G")),
        None
    );
    assert_eq!(
        filter.metadata_decision(&variant("rs2", "2", 20, "C", "T")),
        Some(false)
    );
}
