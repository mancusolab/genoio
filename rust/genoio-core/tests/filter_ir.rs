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
        qual: None,
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

#[test]
fn partial_filter_decision_accepts_metadata_true_or_without_genotypes() {
    let filter = genoio_core::VariantFilter::from_json_value(serde_json::json!({
        "op": "or",
        "left": {"op": "predicate", "name": "qual", "params": {"min": 20.0}},
        "right": {"op": "predicate", "name": "maf", "params": {"max": 0.05}}
    }))
    .expect("filter IR should deserialize");

    let mut high_qual = variant("high_qual", "1", 10, "A", "G");
    high_qual.qual = Some(30.0);

    assert_eq!(
        filter.partial_decision(&high_qual),
        genoio_core::PartialFilterDecision::Accept
    );
}

#[test]
fn partial_filter_decision_rejects_metadata_false_and_without_genotypes() {
    let filter = genoio_core::VariantFilter::from_json_value(serde_json::json!({
        "op": "and",
        "left": {"op": "predicate", "name": "qual", "params": {"min": 20.0}},
        "right": {"op": "predicate", "name": "maf", "params": {"max": 0.05}}
    }))
    .expect("filter IR should deserialize");

    let mut low_qual = variant("low_qual", "1", 10, "A", "G");
    low_qual.qual = Some(10.0);

    assert_eq!(
        filter.partial_decision(&low_qual),
        genoio_core::PartialFilterDecision::Reject
    );
}

#[test]
fn partial_filter_decision_requests_genotypes_when_metadata_cannot_decide() {
    let filter = genoio_core::VariantFilter::from_json_value(serde_json::json!({
        "op": "or",
        "left": {"op": "predicate", "name": "qual", "params": {"min": 20.0}},
        "right": {"op": "predicate", "name": "maf", "params": {"max": 0.05}}
    }))
    .expect("filter IR should deserialize");

    let mut low_qual = variant("low_qual", "1", 10, "A", "G");
    low_qual.qual = Some(10.0);

    assert_eq!(
        filter.partial_decision(&low_qual),
        genoio_core::PartialFilterDecision::NeedGenotypes
    );
}

#[test]
fn qual_predicate_uses_metadata_drop_decisions() {
    let filter = genoio_core::VariantFilter::from_json_value(serde_json::json!({
        "op": "predicate",
        "name": "qual",
        "params": {"min": 20.0, "max": 40.0}
    }))
    .expect("filter IR should deserialize");

    let mut low = variant("low", "1", 10, "A", "G");
    low.qual = Some(10.0);
    assert_eq!(filter.metadata_decision(&low), Some(false));
    assert!(!filter.evaluate(&low, None));

    let mut retained = variant("retained", "1", 20, "A", "G");
    retained.qual = Some(30.0);
    assert_eq!(filter.metadata_decision(&retained), Some(true));
    assert!(filter.evaluate(&retained, None));

    let missing = variant("missing", "1", 30, "A", "G");
    assert_eq!(filter.metadata_decision(&missing), Some(false));
    assert!(!filter.evaluate(&missing, None));
}

#[test]
fn concrete_region_pushdown_is_extracted_only_from_safe_expression_shapes() {
    let and_filter = genoio_core::VariantFilter::from_json_value(serde_json::json!({
        "op": "and",
        "left": {"op": "predicate", "name": "region", "params": {"value": "1:10-20"}},
        "right": {"op": "predicate", "name": "snp", "params": {}}
    }))
    .expect("filter IR should deserialize");

    assert_eq!(
        and_filter.concrete_region_pushdown(),
        Some(genoio_core::RegionPredicate {
            chrom: "1".to_string(),
            start: 10,
            end: 20,
        })
    );

    let or_filter = genoio_core::VariantFilter::from_json_value(serde_json::json!({
        "op": "or",
        "left": {"op": "predicate", "name": "region", "params": {"value": "1:10-20"}},
        "right": {"op": "predicate", "name": "chrom", "params": {"value": "2"}}
    }))
    .expect("filter IR should deserialize");
    assert_eq!(or_filter.concrete_region_pushdown(), None);

    let not_filter = genoio_core::VariantFilter::from_json_value(serde_json::json!({
        "op": "not",
        "expr": {"op": "predicate", "name": "region", "params": {"value": "1:10-20"}}
    }))
    .expect("filter IR should deserialize");
    assert_eq!(not_filter.concrete_region_pushdown(), None);
}

#[test]
fn genotype_stats_preserve_integer_mac_beyond_f32_exact_range() {
    let n_called = 16_777_217_usize;
    let values = vec![1.0_f32; n_called];
    let missing = vec![false; n_called];

    let stats =
        genoio_core::compute_variant_stats(&values, &missing).expect("stats should compute");

    assert_eq!(stats.n_called, u32::try_from(n_called).unwrap());
    assert_eq!(stats.mac, Some(u32::try_from(n_called).unwrap()));
    assert_eq!(stats.af, Some(0.5));
    assert_eq!(stats.maf, Some(0.5));
}
