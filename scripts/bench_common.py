# pattern: Mixed (unavoidable)
# Reason: Benchmark timing, output validation, and iterator cleanup share one support boundary.

from __future__ import annotations

import argparse
import statistics
import time
from collections.abc import Callable
from typing import Any

import numpy as np
from scipy import sparse as scipy_sparse


def read_first_block(dataset: Any, size: int, **read_options: object) -> Any:
    """Return one block and close its persistent reader before returning."""
    with dataset.iter_blocks(size, **read_options) as blocks:
        return next(blocks)


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
    if scipy_sparse.issparse(matrix):
        data = np.asarray(matrix.data)
        return {
            "shape": tuple(matrix.shape),
            "dtype": str(data.dtype),
            "sum": float(data.sum()) if data.size else 0.0,
            "missing": int(np.isnan(data).sum()) if np.issubdtype(data.dtype, np.floating) else 0,
        }
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
