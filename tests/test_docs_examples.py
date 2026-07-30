# pattern: Functional Core

import re
from pathlib import Path


def test_docs_do_not_use_removed_iter_blocks_block_size_keyword():
    docs_dir = Path("docs")
    offenders = [path for path in docs_dir.rglob("*.md") if "block_size=" in path.read_text()]

    assert offenders == []


def test_iter_blocks_docs_use_native_context_manager_for_early_exit():
    api_source = Path("src/genoio/_api.py").read_text()

    assert "contextlib.closing" not in api_source
    assert "with dataset.iter_blocks" in api_source


def test_iter_blocks_docs_explain_generator_compatibility_change():
    reading_docs = Path("docs/api/reading.md").read_text()

    assert "`send()` and `throw()`" in reading_docs
    assert "There is no direct `throw()` equivalent." in reading_docs


def test_public_iter_blocks_examples_use_context_manager():
    public_docs = (
        Path("README.md"),
        Path("docs/index.md"),
        Path("docs/examples/gwas.md"),
        Path("docs/filtering.md"),
        Path("docs/faq.md"),
    )
    direct_loops = re.compile(r"^\s*for .+ in .+\.iter_blocks\(", re.MULTILINE)

    offenders = [path for path in public_docs if direct_loops.search(path.read_text())]

    assert offenders == []
