# pattern: Imperative Shell

from __future__ import annotations

from pathlib import Path


def test_performance_docs_record_phase_6_packed_batch_decision_and_provenance() -> None:
    performance_doc = Path("docs/performance.md").read_text()
    normalized_doc = " ".join(performance_doc.split())

    assert "Phase 6 PLINK2 packed-batch decision run" in normalized_doc
    assert "3bb767085c43c8a39687fa93e4b238c305d3c5bc" in normalized_doc
    assert "0.0100 s for `genoio` and 0.0069 s for `pgenlib`" in normalized_doc
    assert "0.0941 s for `genoio` and 0.0497 s for `pgenlib`" in normalized_doc
    assert "keep packed batches for unfiltered dense source windows" in normalized_doc
