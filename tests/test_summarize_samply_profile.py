# pattern: Mixed

from __future__ import annotations

import gzip
import importlib.util
import json
from pathlib import Path
from types import ModuleType


def load_script(module_name: str) -> ModuleType:
    path = Path(__file__).resolve().parents[1] / "scripts" / f"{module_name}.py"
    spec = importlib.util.spec_from_file_location(module_name, path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


summarize_samply_profile = load_script("summarize_samply_profile")


def synthetic_profile() -> dict[str, object]:
    return {
        "libs": [
            {
                "name": "libtest.dylib",
                "debugName": "libtest.dylib",
                "codeId": "ABCDEF",
            }
        ],
        "threads": [
            {
                "name": "python",
                "processName": "python",
                "stringArray": ["libtest.dylib", "caller", "0x1010", "plain_leaf"],
                "resourceTable": {
                    "length": 1,
                    "lib": [0],
                    "name": [0],
                    "host": [None],
                    "type": [1],
                },
                "funcTable": {
                    "length": 3,
                    "name": [1, 2, 3],
                    "resource": [0, 0, 0],
                    "isJS": [False, False, False],
                    "relevantForJS": [False, False, False],
                    "fileName": [None, None, None],
                    "lineNumber": [None, None, None],
                    "columnNumber": [None, None, None],
                },
                "frameTable": {
                    "length": 3,
                    "func": [0, 1, 2],
                    "nativeSymbol": [None, None, None],
                    "address": [None, None, None],
                    "category": [0, 0, 0],
                    "subcategory": [0, 0, 0],
                    "line": [None, None, None],
                    "column": [None, None, None],
                    "innerWindowID": [None, None, None],
                    "inlineDepth": [0, 0, 0],
                },
                "stackTable": {
                    "length": 3,
                    "prefix": [None, 0, 0],
                    "frame": [0, 1, 2],
                },
                "samples": {
                    "length": 3,
                    "stack": [1, 1, 2],
                    "weight": [1, 2, 1],
                },
            }
        ],
    }


def synthetic_symbols() -> dict[str, object]:
    return {
        "string_table": ["resolved_leaf"],
        "data": [
            {
                "debug_name": "libtest.dylib",
                "code_id": "ABCDEF",
                "symbol_table": [
                    {
                        "rva": 0x1000,
                        "size": 0x100,
                        "symbol": 0,
                    }
                ],
            }
        ],
    }


def test_summarize_profile_resolves_sidecar_symbols() -> None:
    summary = summarize_samply_profile.summarize_profile(
        synthetic_profile(),
        synthetic_symbols(),
        limit=10,
        pattern=None,
        thread_name="python",
    )

    assert summary.thread_name == "python"
    assert summary.total_weight == 4
    assert [(row.name, row.count) for row in summary.inclusive_rows] == [
        ("caller", 4),
        ("resolved_leaf", 3),
        ("plain_leaf", 1),
    ]
    assert [(row.name, row.count) for row in summary.self_rows] == [
        ("resolved_leaf", 3),
        ("plain_leaf", 1),
    ]


def test_summarize_profile_filters_pattern() -> None:
    summary = summarize_samply_profile.summarize_profile(
        synthetic_profile(),
        synthetic_symbols(),
        limit=10,
        pattern="leaf",
        thread_name="python",
    )

    assert [row.name for row in summary.inclusive_rows] == ["resolved_leaf", "plain_leaf"]
    assert [row.name for row in summary.self_rows] == ["resolved_leaf", "plain_leaf"]


def test_load_profile_and_matching_sidecar(tmp_path: Path) -> None:
    profile_path = tmp_path / "sample.profile.json.gz"
    sidecar_path = tmp_path / "sample.profile.json.syms.json"
    with gzip.open(profile_path, "wt", encoding="utf-8") as handle:
        json.dump(synthetic_profile(), handle)
    sidecar_path.write_text(json.dumps(synthetic_symbols()), encoding="utf-8")

    profile = summarize_samply_profile.load_profile(profile_path)
    symbols = summarize_samply_profile.load_matching_symbols(profile_path)

    assert profile["threads"][0]["name"] == "python"
    assert symbols["data"][0]["debug_name"] == "libtest.dylib"


def test_load_matching_symbols_returns_none_without_sidecar(tmp_path: Path) -> None:
    profile_path = tmp_path / "sample.profile.json.gz"
    with gzip.open(profile_path, "wt", encoding="utf-8") as handle:
        json.dump(synthetic_profile(), handle)

    assert summarize_samply_profile.load_matching_symbols(profile_path) is None
