# pattern: Imperative Shell

from __future__ import annotations

from pathlib import Path


def test_performance_docs_record_direct_fill_decision_and_provenance() -> None:
    performance_doc = Path("docs/performance.md").read_text()
    normalized_doc = " ".join(performance_doc.split())

    assert "direct sample-major fill" in normalized_doc
    assert "0.0145-0.0147 s" in normalized_doc
    assert "135f3d2a165882263ed3520872998473bfd9615b" in normalized_doc
    assert "removing the variant-major accumulation and transpose copy" in normalized_doc
