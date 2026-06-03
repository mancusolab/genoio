from __future__ import annotations

import argparse
import statistics
import time
from collections.abc import Callable, Iterator
from contextlib import contextmanager
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any

import numpy as np


def positive_int(value: str) -> int:
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("value must be a positive integer")
    return parsed


def nonnegative_float(value: str) -> float:
    parsed = float(value)
    if parsed < 0:
        raise argparse.ArgumentTypeError("value must be nonnegative")
    return parsed


def matrix_summary(matrix: Any) -> dict[str, Any]:
    array = np.asarray(matrix)
    return {
        "shape": tuple(array.shape),
        "dtype": str(array.dtype),
        "sum": float(np.nansum(array)),
        "missing": int(np.isnan(array).sum()) if np.issubdtype(array.dtype, np.floating) else 0,
    }


def benchmark(name: str, fn: Callable[[], Any], repeats: int) -> Any:
    first = fn()
    summary = matrix_summary(first)
    times = []
    for _ in range(repeats):
        start = time.perf_counter()
        observed = fn()
        times.append(time.perf_counter() - start)
        observed_summary = matrix_summary(observed)
        if observed_summary["shape"] != summary["shape"]:
            raise RuntimeError(f"{name} shape changed between repeats")
    print_result(name, summary, times)
    return first


def print_result(name: str, summary: dict[str, Any], times: list[float]) -> None:
    print(name)
    print(
        "  matrix",
        f"shape={summary['shape']}",
        f"dtype={summary['dtype']}",
        f"sum={summary['sum']:.6g}",
        f"missing={summary['missing']}",
    )
    print(
        "  time",
        f"median={statistics.median(times):.4f}s",
        f"min={min(times):.4f}s",
        "runs=" + " ".join(f"{value:.4f}" for value in times),
    )


def compare_summaries(left_name: str, left: Any, right_name: str, right: Any) -> None:
    left_array = np.asarray(left)
    right_array = np.asarray(right)
    print("comparison")
    print(f"  {left_name}.shape={left_array.shape} {right_name}.shape={right_array.shape}")
    if left_array.shape != right_array.shape:
        print("  skipped value comparison: shapes differ")
        return
    equal = np.allclose(left_array, right_array, equal_nan=True)
    max_abs_diff = float(np.nanmax(np.abs(left_array - right_array))) if left_array.size else 0.0
    print(f"  allclose={equal} max_abs_diff={max_abs_diff:.6g}")


@contextmanager
def plink2_prefix_with_uncompressed_pvar(prefix: Path) -> Iterator[Path]:
    if prefix.with_suffix(".pvar").exists():
        yield prefix
        return
    pvar_zst = prefix.with_suffix(".pvar.zst")
    if not pvar_zst.exists():
        yield prefix
        return

    import shutil
    import subprocess

    zstd = shutil.which("zstd")
    if zstd is None:
        raise RuntimeError(f"{prefix}.pvar is missing and zstd is not available to decompress {pvar_zst}")

    with TemporaryDirectory(prefix="genoio-plink2-bench-") as tmpdir:
        tmp_prefix = Path(tmpdir) / prefix.name
        for suffix in (".pgen", ".psam"):
            source = prefix.with_suffix(suffix).resolve()
            target = tmp_prefix.with_suffix(suffix)
            try:
                target.symlink_to(source)
                continue
            except OSError:
                pass
            shutil.copy2(source, target)
        with tmp_prefix.with_suffix(".pvar").open("wb") as out:
            subprocess.run([zstd, "-dc", str(pvar_zst)], stdout=out, check=True)
        yield tmp_prefix
