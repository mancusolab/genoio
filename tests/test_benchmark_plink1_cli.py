# pattern: Mixed

from __future__ import annotations

import sys
from pathlib import Path
from types import ModuleType
from typing import Any, cast

import numpy as np
from script_loader import load_benchmark_script

benchmark_plink1 = cast(Any, load_benchmark_script("benchmark_plink1"))


class _FakePlinkArray:
    values = np.array([[0.0, 1.0]], dtype=np.float32)

    def isel(self, *, variant: slice) -> _FakePlinkArray:
        assert variant == slice(0, 7)
        return self


def test_parse_args_accepts_pandas_plink_options(monkeypatch) -> None:
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "benchmark_plink1.py",
            "--prefix",
            "data/example",
            "--backend",
            "pandas_plink",
            "--max-variants",
            "7",
            "--repeats",
            "2",
            "--pandas-ref",
            "a1",
            "--no-compare",
        ],
    )

    args = benchmark_plink1.parse_args()

    assert args.prefix == Path("data/example")
    assert args.backend == "pandas_plink"
    assert args.max_variants == 7
    assert args.repeats == 2
    assert args.pandas_ref == "a1"
    assert args.no_compare is True


def test_read_pandas_plink_uses_string_paths_and_patches_dask(monkeypatch) -> None:
    observed: dict[str, Any] = {"patched": False}
    fake_module = ModuleType("pandas_plink")

    def read_plink1_bin(bed: str, bim: str, fam: str, *, verbose: bool, ref: str) -> _FakePlinkArray:
        observed["args"] = (bed, bim, fam, verbose, ref)
        return _FakePlinkArray()

    fake_module.read_plink1_bin = read_plink1_bin  # ty: ignore[unresolved-attribute]
    monkeypatch.setitem(sys.modules, "pandas_plink", fake_module)
    monkeypatch.setattr(benchmark_plink1, "_patch_dask_memmap_tokenization", lambda: observed.update(patched=True))
    args = benchmark_plink1.argparse.Namespace(prefix=Path("data/example"), max_variants=7, pandas_ref="a0")

    matrix = benchmark_plink1.read_pandas_plink(args)

    assert observed["patched"] is True
    assert observed["args"] == (
        "data/example.bed",
        "data/example.bim",
        "data/example.fam",
        False,
        "a0",
    )
    np.testing.assert_array_equal(matrix, np.array([[0.0, 1.0]], dtype=np.float32))
