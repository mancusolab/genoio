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

#[derive(Clone)]
enum LeafPredicate {
    Chrom1,
    IdRs1,
    QualMin20,
    MafMin40,
}

#[derive(Clone)]
enum RefExpr {
    Leaf(LeafPredicate),
    And(Box<RefExpr>, Box<RefExpr>),
    Or(Box<RefExpr>, Box<RefExpr>),
    Not(Box<RefExpr>),
}

#[derive(Clone, Copy)]
struct LeafValues {
    chrom1: bool,
    id_rs1: bool,
    qual_min20: bool,
    maf_min40: bool,
}

fn generated_boolean_expressions() -> Vec<RefExpr> {
    let leaves = vec![
        RefExpr::Leaf(LeafPredicate::Chrom1),
        RefExpr::Leaf(LeafPredicate::IdRs1),
        RefExpr::Leaf(LeafPredicate::QualMin20),
        RefExpr::Leaf(LeafPredicate::MafMin40),
    ];
    let mut expressions = leaves.clone();
    expressions.extend(
        leaves
            .iter()
            .cloned()
            .map(|expr| RefExpr::Not(Box::new(expr))),
    );
    for left in &leaves {
        for right in &leaves {
            expressions.push(RefExpr::And(
                Box::new(left.clone()),
                Box::new(right.clone()),
            ));
            expressions.push(RefExpr::Or(Box::new(left.clone()), Box::new(right.clone())));
        }
    }

    let one_level = expressions.clone();
    for left in &one_level {
        for right in &one_level {
            expressions.push(RefExpr::And(
                Box::new(left.clone()),
                Box::new(right.clone()),
            ));
            expressions.push(RefExpr::Or(Box::new(left.clone()), Box::new(right.clone())));
        }
    }
    expressions.extend(
        one_level
            .into_iter()
            .map(|expr| RefExpr::Not(Box::new(expr))),
    );
    expressions
}

fn ref_expr_to_json(expr: &RefExpr) -> serde_json::Value {
    match expr {
        RefExpr::Leaf(LeafPredicate::Chrom1) => serde_json::json!({
            "op": "predicate",
            "name": "chrom",
            "params": {"value": "1"}
        }),
        RefExpr::Leaf(LeafPredicate::IdRs1) => serde_json::json!({
            "op": "predicate",
            "name": "id_in",
            "params": {"values": ["rs1"]}
        }),
        RefExpr::Leaf(LeafPredicate::QualMin20) => serde_json::json!({
            "op": "predicate",
            "name": "qual",
            "params": {"min": 20.0}
        }),
        RefExpr::Leaf(LeafPredicate::MafMin40) => serde_json::json!({
            "op": "predicate",
            "name": "maf",
            "params": {"min": 0.4}
        }),
        RefExpr::And(left, right) => serde_json::json!({
            "op": "and",
            "left": ref_expr_to_json(left),
            "right": ref_expr_to_json(right)
        }),
        RefExpr::Or(left, right) => serde_json::json!({
            "op": "or",
            "left": ref_expr_to_json(left),
            "right": ref_expr_to_json(right)
        }),
        RefExpr::Not(expr) => serde_json::json!({
            "op": "not",
            "expr": ref_expr_to_json(expr)
        }),
    }
}

fn eval_ref_expr(expr: &RefExpr, values: LeafValues) -> bool {
    match expr {
        RefExpr::Leaf(LeafPredicate::Chrom1) => values.chrom1,
        RefExpr::Leaf(LeafPredicate::IdRs1) => values.id_rs1,
        RefExpr::Leaf(LeafPredicate::QualMin20) => values.qual_min20,
        RefExpr::Leaf(LeafPredicate::MafMin40) => values.maf_min40,
        RefExpr::And(left, right) => eval_ref_expr(left, values) && eval_ref_expr(right, values),
        RefExpr::Or(left, right) => eval_ref_expr(left, values) || eval_ref_expr(right, values),
        RefExpr::Not(expr) => !eval_ref_expr(expr, values),
    }
}

fn eval_ref_metadata_decision(expr: &RefExpr, values: LeafValues) -> Option<bool> {
    match expr {
        RefExpr::Leaf(LeafPredicate::Chrom1) => Some(values.chrom1),
        RefExpr::Leaf(LeafPredicate::IdRs1) => Some(values.id_rs1),
        RefExpr::Leaf(LeafPredicate::QualMin20) => Some(values.qual_min20),
        RefExpr::Leaf(LeafPredicate::MafMin40) => None,
        RefExpr::And(left, right) => match (
            eval_ref_metadata_decision(left, values),
            eval_ref_metadata_decision(right, values),
        ) {
            (Some(false), _) | (_, Some(false)) => Some(false),
            (Some(true), Some(true)) => Some(true),
            _ => None,
        },
        RefExpr::Or(left, right) => match (
            eval_ref_metadata_decision(left, values),
            eval_ref_metadata_decision(right, values),
        ) {
            (Some(true), _) | (_, Some(true)) => Some(true),
            (Some(false), Some(false)) => Some(false),
            _ => None,
        },
        RefExpr::Not(expr) => eval_ref_metadata_decision(expr, values).map(|decision| !decision),
    }
}

fn variant_for_values(values: LeafValues) -> genoio_core::VariantRecord {
    let mut record = variant(
        if values.id_rs1 { "rs1" } else { "rs2" },
        if values.chrom1 { "1" } else { "2" },
        10,
        "A",
        "G",
    );
    record.qual = Some(if values.qual_min20 { 30.0 } else { 10.0 });
    record
}

fn stats_for_values(values: LeafValues) -> genoio_core::VariantStats {
    genoio_core::VariantStats {
        af: Some(if values.maf_min40 { 0.5 } else { 0.1 }),
        maf: Some(if values.maf_min40 { 0.5 } else { 0.1 }),
        mac: Some(if values.maf_min40 { 4.0 } else { 1.0 }),
        missing_rate: 0.0,
        n_called: 4,
        polymorphic: true,
    }
}

fn all_leaf_values() -> Vec<LeafValues> {
    let mut values = Vec::new();
    for bits in 0_u8..16 {
        values.push(LeafValues {
            chrom1: bits & 0b0001 != 0,
            id_rs1: bits & 0b0010 != 0,
            qual_min20: bits & 0b0100 != 0,
            maf_min40: bits & 0b1000 != 0,
        });
    }
    values
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
fn boolean_filter_ir_matches_reference_truth_table_for_nested_expressions() {
    for expr in generated_boolean_expressions() {
        let filter = genoio_core::VariantFilter::from_json_value(ref_expr_to_json(&expr))
            .expect("generated filter should parse");

        for values in all_leaf_values() {
            let variant = variant_for_values(values);
            let stats = stats_for_values(values);

            assert_eq!(
                filter.evaluate(&variant, Some(&stats)),
                eval_ref_expr(&expr, values)
            );
        }
    }
}

#[test]
fn partial_filter_decisions_match_reference_truth_table_when_metadata_is_decisive() {
    for expr in generated_boolean_expressions() {
        let filter = genoio_core::VariantFilter::from_json_value(ref_expr_to_json(&expr))
            .expect("generated filter should parse");

        for values in all_leaf_values() {
            let variant = variant_for_values(values);
            let expected = eval_ref_metadata_decision(&expr, values);
            assert_eq!(filter.metadata_decision(&variant), expected);
            assert_eq!(
                filter.partial_decision(&variant),
                match expected {
                    Some(true) => genoio_core::PartialFilterDecision::Accept,
                    Some(false) => genoio_core::PartialFilterDecision::Reject,
                    None => genoio_core::PartialFilterDecision::NeedGenotypes,
                }
            );
        }
    }
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
    assert!(!filter.is_genotype_stats_only());
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
fn genotype_stats_only_filters_evaluate_without_variant_metadata() {
    let filter = genoio_core::VariantFilter::from_json_value(serde_json::json!({
        "op": "and",
        "left": {"op": "predicate", "name": "maf", "params": {"min": 0.2}},
        "right": {"op": "predicate", "name": "missing_rate", "params": {"max": 0.1}}
    }))
    .expect("filter IR should deserialize");
    let retained = genoio_core::VariantStats {
        af: Some(0.25),
        maf: Some(0.25),
        mac: Some(4.0),
        missing_rate: 0.0,
        n_called: 8,
        polymorphic: true,
    };
    let rejected = genoio_core::VariantStats {
        maf: Some(0.1),
        ..retained
    };

    assert!(filter.requires_genotype_stats());
    assert!(filter.is_genotype_stats_only());
    assert_eq!(filter.evaluate_genotype_stats(&retained), Some(true));
    assert_eq!(filter.evaluate_genotype_stats(&rejected), Some(false));
}

#[test]
fn genotype_filter_plan_matches_genotype_stats_only_evaluation() {
    let filter = genoio_core::VariantFilter::from_json_value(serde_json::json!({
        "op": "and",
        "left": {"op": "predicate", "name": "maf", "params": {"min": 0.2, "max": 0.4}},
        "right": {"op": "predicate", "name": "missing_rate", "params": {"max": 0.1}}
    }))
    .expect("filter IR should deserialize");
    let plan = filter.genotype_filter_plan();

    for stats in [
        genoio_core::VariantStats {
            af: Some(0.25),
            maf: Some(0.25),
            mac: Some(4.0),
            missing_rate: 0.0,
            n_called: 8,
            polymorphic: true,
        },
        genoio_core::VariantStats {
            af: Some(0.45),
            maf: Some(0.45),
            mac: Some(7.0),
            missing_rate: 0.0,
            n_called: 8,
            polymorphic: true,
        },
        genoio_core::VariantStats {
            af: Some(0.25),
            maf: Some(0.25),
            mac: Some(4.0),
            missing_rate: 0.2,
            n_called: 8,
            polymorphic: true,
        },
    ] {
        assert_eq!(
            plan.evaluate_stats(&stats),
            filter.evaluate_genotype_stats(&stats)
        );
    }
}

#[test]
fn polymorphic_filters_compile_to_specialized_genotype_plan() {
    let plain = genoio_core::VariantFilter::from_json_value(serde_json::json!({
        "op": "predicate",
        "name": "polymorphic",
        "params": {}
    }))
    .expect("filter IR should deserialize");
    let metadata_gated = genoio_core::VariantFilter::from_json_value(serde_json::json!({
        "op": "and",
        "left": {"op": "predicate", "name": "biallelic", "params": {}},
        "right": {"op": "predicate", "name": "polymorphic", "params": {}}
    }))
    .expect("filter IR should deserialize");
    let conjoined = genoio_core::VariantFilter::from_json_value(serde_json::json!({
        "op": "and",
        "left": {"op": "predicate", "name": "polymorphic", "params": {}},
        "right": {"op": "predicate", "name": "missing_rate", "params": {"max": 0.1}}
    }))
    .expect("filter IR should deserialize");

    assert_eq!(
        plain.genotype_filter_plan(),
        genoio_core::GenotypeFilterPlan::Polymorphic
    );
    assert_eq!(
        metadata_gated.genotype_filter_plan(),
        genoio_core::GenotypeFilterPlan::Polymorphic
    );
    assert_eq!(
        conjoined.genotype_filter_plan(),
        genoio_core::GenotypeFilterPlan::Conjunction(genoio_core::GenotypeFilterConjunction {
            polymorphic: true,
            mac_min: None,
            mac_max: None,
            maf_min: None,
            maf_max: None,
            missing_rate_max: Some(0.1),
        })
    );
}

#[test]
fn range_filters_compile_to_specialized_genotype_plans() {
    let stats = genoio_core::VariantStats {
        af: Some(0.125),
        maf: Some(0.125),
        mac: Some(4.0),
        missing_rate: 0.02,
        n_called: 16,
        polymorphic: true,
    };

    let mac_filter = genoio_core::VariantFilter::from_json_value(serde_json::json!({
        "op": "predicate",
        "name": "mac",
        "params": {"min": 2, "max": 8}
    }))
    .expect("mac filter should parse");
    let mac_plan = mac_filter.genotype_filter_plan();
    assert_eq!(
        mac_plan,
        genoio_core::GenotypeFilterPlan::MacRange {
            min: Some(2),
            max: Some(8)
        }
    );
    assert_eq!(mac_plan.evaluate_stats(&stats), Some(true));

    let maf_filter = genoio_core::VariantFilter::from_json_value(serde_json::json!({
        "op": "and",
        "left": {"op": "predicate", "name": "chrom", "params": {"value": "1"}},
        "right": {"op": "predicate", "name": "maf", "params": {"min": 0.1, "max": 0.2}}
    }))
    .expect("mixed metadata/maf filter should parse");
    let maf_plan = maf_filter.genotype_filter_plan();
    assert_eq!(
        maf_plan,
        genoio_core::GenotypeFilterPlan::MafRange {
            min: Some(0.1),
            max: Some(0.2)
        }
    );
    assert_eq!(maf_plan.evaluate_stats(&stats), Some(true));

    let missing_filter = genoio_core::VariantFilter::from_json_value(serde_json::json!({
        "op": "predicate",
        "name": "missing_rate",
        "params": {"max": 0.01}
    }))
    .expect("missing-rate filter should parse");
    let missing_plan = missing_filter.genotype_filter_plan();
    assert_eq!(
        missing_plan,
        genoio_core::GenotypeFilterPlan::MissingRateMax { max: 0.01 }
    );
    assert_eq!(missing_plan.evaluate_stats(&stats), Some(false));
}

#[test]
fn conjunction_filters_compile_to_specialized_genotype_plan() {
    let filter = genoio_core::VariantFilter::from_json_value(serde_json::json!({
        "op": "and",
        "left": {"op": "predicate", "name": "mac", "params": {"min": 2}},
        "right": {"op": "predicate", "name": "missing_rate", "params": {"max": 0.05}}
    }))
    .expect("conjoined genotype filter should parse");

    let plan = filter.genotype_filter_plan();
    assert_eq!(
        plan,
        genoio_core::GenotypeFilterPlan::Conjunction(genoio_core::GenotypeFilterConjunction {
            polymorphic: false,
            mac_min: Some(2),
            mac_max: None,
            maf_min: None,
            maf_max: None,
            missing_rate_max: Some(0.05),
        })
    );
    assert_eq!(
        plan.evaluate_stats(&genoio_core::VariantStats {
            af: Some(0.125),
            maf: Some(0.125),
            mac: Some(4.0),
            missing_rate: 0.02,
            n_called: 16,
            polymorphic: true,
        }),
        Some(true)
    );
    assert_eq!(
        plan.evaluate_stats(&genoio_core::VariantStats {
            af: Some(0.125),
            maf: Some(0.125),
            mac: Some(4.0),
            missing_rate: 0.10,
            n_called: 16,
            polymorphic: true,
        }),
        Some(false)
    );
}

#[test]
fn unsupported_boolean_filters_keep_generic_genotype_plan() {
    let filter = genoio_core::VariantFilter::from_json_value(serde_json::json!({
        "op": "or",
        "left": {"op": "predicate", "name": "mac", "params": {"min": 2}},
        "right": {"op": "predicate", "name": "missing_rate", "params": {"max": 0.05}}
    }))
    .expect("or filter should parse");

    let plan = filter.genotype_filter_plan();
    assert_eq!(plan, genoio_core::GenotypeFilterPlan::Generic);
    assert_eq!(
        plan.evaluate_stats(&genoio_core::VariantStats {
            af: Some(0.125),
            maf: Some(0.125),
            mac: Some(4.0),
            missing_rate: 0.02,
            n_called: 16,
            polymorphic: true,
        }),
        None
    );
}

#[test]
fn dosage_values_detect_polymorphic_without_full_variant_stats() {
    assert!(
        genoio_core::is_dosage_polymorphic(&[0.0, 1.0, 0.0], &[false, false, false])
            .expect("valid dosages should evaluate")
    );
    assert!(
        !genoio_core::is_dosage_polymorphic(&[0.0, 0.0, 0.0], &[false, false, true])
            .expect("valid dosages should evaluate")
    );
    assert!(
        !genoio_core::is_dosage_polymorphic(&[2.0, 2.0], &[false, false])
            .expect("valid dosages should evaluate")
    );
    assert!(!genoio_core::is_dosage_polymorphic(&[1.0], &[true])
        .expect("missing dosages should not count"));
    assert!(
        genoio_core::is_dosage_polymorphic(&[0.1, 1.9], &[false, false])
            .expect("fractional dosages should evaluate")
    );
    assert!(genoio_core::is_dosage_polymorphic(&[3.0], &[false]).is_err());
}

#[test]
fn genotype_stats_only_detection_handles_boolean_composition() {
    let cases = [
        (
            serde_json::json!({
                "op": "or",
                "left": {"op": "predicate", "name": "maf", "params": {"min": 0.2}},
                "right": {"op": "predicate", "name": "missing_rate", "params": {"max": 0.1}}
            }),
            true,
        ),
        (
            serde_json::json!({
                "op": "not",
                "expr": {"op": "predicate", "name": "polymorphic", "params": {}}
            }),
            true,
        ),
        (
            serde_json::json!({
                "op": "or",
                "left": {"op": "predicate", "name": "id_in", "params": {"values": ["rs1"]}},
                "right": {"op": "predicate", "name": "maf", "params": {"min": 0.2}}
            }),
            false,
        ),
    ];

    for (expr, expected) in cases {
        let filter = genoio_core::VariantFilter::from_json_value(expr)
            .expect("filter IR should deserialize");

        assert_eq!(filter.is_genotype_stats_only(), expected);
    }
}

#[test]
fn mixed_metadata_filters_do_not_evaluate_from_genotype_stats_only() {
    let filter = genoio_core::VariantFilter::from_json_value(serde_json::json!({
        "op": "and",
        "left": {"op": "predicate", "name": "snp", "params": {}},
        "right": {"op": "predicate", "name": "maf", "params": {"min": 0.2}}
    }))
    .expect("filter IR should deserialize");
    let stats = genoio_core::VariantStats {
        af: Some(0.25),
        maf: Some(0.25),
        mac: Some(4.0),
        missing_rate: 0.0,
        n_called: 8,
        polymorphic: true,
    };

    assert!(filter.requires_genotype_stats());
    assert!(!filter.is_genotype_stats_only());
    assert_eq!(filter.evaluate_genotype_stats(&stats), None);
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
fn concrete_region_pushdown_intersects_conjoined_regions() {
    let overlapping = genoio_core::VariantFilter::from_json_value(serde_json::json!({
        "op": "and",
        "left": {"op": "predicate", "name": "region", "params": {"value": "1:10-30"}},
        "right": {"op": "predicate", "name": "region", "params": {"value": "1:20-40"}}
    }))
    .expect("filter IR should deserialize");

    assert_eq!(
        overlapping.concrete_region_pushdown(),
        Some(genoio_core::RegionPredicate {
            chrom: "1".to_string(),
            start: 20,
            end: 30,
        })
    );

    let disjoint = genoio_core::VariantFilter::from_json_value(serde_json::json!({
        "op": "and",
        "left": {"op": "predicate", "name": "region", "params": {"value": "1:10-20"}},
        "right": {"op": "predicate", "name": "region", "params": {"value": "1:30-40"}}
    }))
    .expect("filter IR should deserialize");

    assert_eq!(disjoint.concrete_region_pushdown(), None);
    assert!(disjoint.is_always_false());
}

#[test]
fn concrete_region_pushdown_combines_chrom_with_region() {
    let matching_chrom = genoio_core::VariantFilter::from_json_value(serde_json::json!({
        "op": "and",
        "left": {"op": "predicate", "name": "chrom", "params": {"value": "1"}},
        "right": {"op": "predicate", "name": "region", "params": {"value": "1:10-20"}}
    }))
    .expect("filter IR should deserialize");

    assert_eq!(
        matching_chrom.concrete_region_pushdown(),
        Some(genoio_core::RegionPredicate {
            chrom: "1".to_string(),
            start: 10,
            end: 20,
        })
    );
    assert!(!matching_chrom.is_always_false());

    let conflicting_chrom = genoio_core::VariantFilter::from_json_value(serde_json::json!({
        "op": "and",
        "left": {"op": "predicate", "name": "chrom", "params": {"value": "2"}},
        "right": {"op": "predicate", "name": "region", "params": {"value": "1:10-20"}}
    }))
    .expect("filter IR should deserialize");

    assert_eq!(conflicting_chrom.concrete_region_pushdown(), None);
    assert!(conflicting_chrom.is_always_false());
}

#[test]
fn conjoined_threshold_predicates_are_tightened_or_rejected() {
    let stats = genoio_core::VariantStats {
        af: Some(0.03),
        maf: Some(0.03),
        mac: Some(3.0),
        missing_rate: 0.02,
        n_called: 100,
        polymorphic: true,
    };
    let variant = variant("rs1", "1", 10, "A", "G");

    let tightened = genoio_core::VariantFilter::from_json_value(serde_json::json!({
        "op": "and",
        "left": {"op": "predicate", "name": "maf", "params": {"min": 0.01}},
        "right": {"op": "predicate", "name": "maf", "params": {"max": 0.05}}
    }))
    .expect("filter IR should deserialize");

    assert!(!tightened.is_always_false());
    assert!(tightened.evaluate(&variant, Some(&stats)));

    let impossible = genoio_core::VariantFilter::from_json_value(serde_json::json!({
        "op": "and",
        "left": {"op": "predicate", "name": "maf", "params": {"min": 0.10}},
        "right": {"op": "predicate", "name": "maf", "params": {"max": 0.05}}
    }))
    .expect("filter IR should deserialize");

    assert!(impossible.is_always_false());
    assert_eq!(
        impossible.partial_decision(&variant),
        genoio_core::PartialFilterDecision::Reject
    );
    assert!(!impossible.evaluate(&variant, Some(&stats)));
}

#[test]
fn id_in_predicates_intersect_and_union() {
    let rs1 = variant("rs1", "1", 10, "A", "G");
    let rs2 = variant("rs2", "1", 20, "A", "G");
    let rs3 = variant("rs3", "1", 30, "A", "G");

    let intersection = genoio_core::VariantFilter::from_json_value(serde_json::json!({
        "op": "and",
        "left": {"op": "predicate", "name": "id_in", "params": {"values": ["rs1", "rs2"]}},
        "right": {"op": "predicate", "name": "id_in", "params": {"values": ["rs2", "rs3"]}}
    }))
    .expect("filter IR should deserialize");

    assert!(!intersection.evaluate(&rs1, None));
    assert!(intersection.evaluate(&rs2, None));
    assert!(!intersection.evaluate(&rs3, None));

    let union = genoio_core::VariantFilter::from_json_value(serde_json::json!({
        "op": "or",
        "left": {"op": "predicate", "name": "id_in", "params": {"values": ["rs1", "rs2"]}},
        "right": {"op": "predicate", "name": "id_in", "params": {"values": ["rs2", "rs3"]}}
    }))
    .expect("filter IR should deserialize");

    assert!(union.evaluate(&rs1, None));
    assert!(union.evaluate(&rs2, None));
    assert!(union.evaluate(&rs3, None));

    let empty_intersection = genoio_core::VariantFilter::from_json_value(serde_json::json!({
        "op": "and",
        "left": {"op": "predicate", "name": "id_in", "params": {"values": ["rs1"]}},
        "right": {"op": "predicate", "name": "id_in", "params": {"values": ["rs2"]}}
    }))
    .expect("filter IR should deserialize");
    assert!(empty_intersection.is_always_false());
}

#[test]
fn genotype_stats_preserve_integer_mac_beyond_f32_exact_range() {
    let n_called = 16_777_217_usize;
    let values = vec![1.0_f32; n_called];
    let missing = vec![false; n_called];

    let stats =
        genoio_core::compute_variant_stats(&values, &missing).expect("stats should compute");

    assert_eq!(stats.n_called, u32::try_from(n_called).unwrap());
    assert_eq!(stats.mac, Some(f64::from(u32::try_from(n_called).unwrap())));
    assert_eq!(stats.af, Some(0.5));
    assert_eq!(stats.maf, Some(0.5));
}
