# pattern: Functional Core

from pathlib import Path


def test_docs_do_not_use_removed_iter_blocks_block_size_keyword():
    docs_dir = Path("docs")
    offenders = [path for path in docs_dir.rglob("*.md") if "block_size=" in path.read_text()]

    assert offenders == []
