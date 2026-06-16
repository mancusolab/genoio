#!/usr/bin/env python
# pattern: Mixed

from __future__ import annotations

import argparse
import bisect
import gzip
import json
import re
from collections import Counter
from pathlib import Path
from typing import Any, NamedTuple


class SummaryRow(NamedTuple):
    name: str
    count: int
    percent: float


class ProfileSummary(NamedTuple):
    thread_name: str
    process_name: str
    total_weight: int
    inclusive_rows: list[SummaryRow]
    self_rows: list[SummaryRow]


class SymbolTable(NamedTuple):
    starts: list[int]
    entries: list[dict[str, Any]]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Summarize top frames from a samply Firefox profile."
    )
    parser.add_argument("profile", type=Path, help="Path to *.profile.json.gz")
    parser.add_argument(
        "--symbols",
        type=Path,
        default=None,
        help="Optional .syms.json sidecar. Defaults to the samply-generated sibling.",
    )
    parser.add_argument("--thread", default="python", help="Thread name to summarize.")
    parser.add_argument("--pattern", default=None, help="Regex for frame names to include.")
    parser.add_argument("--limit", type=positive_int, default=30)
    return parser.parse_args()


def positive_int(value: str) -> int:
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("value must be a positive integer")
    return parsed


def load_profile(path: Path) -> dict[str, Any]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return json.load(handle)


def load_matching_symbols(profile_path: Path) -> dict[str, Any] | None:
    sidecar_path = matching_symbols_path(profile_path)
    if not sidecar_path.exists():
        return None
    with sidecar_path.open(encoding="utf-8") as handle:
        return json.load(handle)


def matching_symbols_path(profile_path: Path) -> Path:
    if profile_path.suffix == ".gz":
        profile_json_path = profile_path.with_suffix("")
    else:
        profile_json_path = profile_path
    return profile_json_path.with_name(f"{profile_json_path.name}.syms.json")


def summarize_profile(
    profile: dict[str, Any],
    symbols: dict[str, Any] | None,
    *,
    limit: int,
    pattern: str | None,
    thread_name: str | None,
) -> ProfileSummary:
    thread = select_thread(profile, thread_name)
    resolver = FrameResolver(profile, thread, symbols)
    name_filter = re.compile(pattern) if pattern is not None else None
    inclusive, self_counts, total_weight = collect_frame_counts(thread, resolver)
    inclusive_rows = top_rows(inclusive, total_weight, limit, name_filter)
    self_rows = top_rows(self_counts, total_weight, limit, name_filter)
    return ProfileSummary(
        thread_name=str(thread.get("name", "")),
        process_name=str(thread.get("processName", "")),
        total_weight=total_weight,
        inclusive_rows=inclusive_rows,
        self_rows=self_rows,
    )


def select_thread(profile: dict[str, Any], thread_name: str | None) -> dict[str, Any]:
    threads = profile.get("threads")
    if not isinstance(threads, list) or not threads:
        raise ValueError("profile contains no threads")
    if thread_name is None:
        return threads[0]
    for thread in threads:
        if thread.get("name") == thread_name:
            return thread
    names = ", ".join(str(thread.get("name", "<unnamed>")) for thread in threads)
    raise ValueError(f"thread {thread_name!r} not found; available threads: {names}")


def collect_frame_counts(
    thread: dict[str, Any],
    resolver: FrameResolver,
) -> tuple[Counter[str], Counter[str], int]:
    samples = thread["samples"]
    stack_table = thread["stackTable"]
    stack_prefix = stack_table["prefix"]
    stack_frame = stack_table["frame"]
    stacks = samples["stack"]
    weights = samples.get("weight")
    if weights is None:
        weights = [1] * int(samples["length"])

    inclusive: Counter[str] = Counter()
    self_counts: Counter[str] = Counter()
    total_weight = 0
    for stack_index, weight in zip(stacks, weights, strict=True):
        if stack_index is None:
            continue
        total_weight += int(weight)
        seen: set[str] = set()
        current = stack_index
        is_leaf = True
        while current is not None:
            frame_name = resolver.frame_name(stack_frame[current])
            if is_leaf:
                self_counts[frame_name] += int(weight)
                is_leaf = False
            if frame_name not in seen:
                inclusive[frame_name] += int(weight)
                seen.add(frame_name)
            current = stack_prefix[current]
    return inclusive, self_counts, total_weight


def top_rows(
    counts: Counter[str],
    total_weight: int,
    limit: int,
    pattern: re.Pattern[str] | None,
) -> list[SummaryRow]:
    rows = []
    for name, count in counts.most_common():
        if pattern is not None and pattern.search(name) is None:
            continue
        percent = 0.0 if total_weight == 0 else 100.0 * count / total_weight
        rows.append(SummaryRow(name=name, count=int(count), percent=percent))
        if len(rows) >= limit:
            break
    return rows


class FrameResolver:
    def __init__(
        self,
        profile: dict[str, Any],
        thread: dict[str, Any],
        symbols: dict[str, Any] | None,
    ) -> None:
        self._libs = profile.get("libs", [])
        self._strings = thread["stringArray"]
        self._resources = thread["resourceTable"]
        self._func_table = thread["funcTable"]
        self._frame_table = thread["frameTable"]
        self._symbols = build_symbol_tables(symbols)

    def frame_name(self, frame_index: int) -> str:
        function_index = self._frame_table["func"][frame_index]
        raw_name = self._strings[self._func_table["name"][function_index]]
        resolved = self._resolve_address(raw_name, function_index)
        return resolved if resolved is not None else raw_name

    def _resolve_address(self, raw_name: str, function_index: int) -> str | None:
        if not raw_name.startswith("0x"):
            return None
        resource_index = self._func_table["resource"][function_index]
        if resource_index is None:
            return None
        lib_index = self._resources["lib"][resource_index]
        if lib_index is None:
            return None
        try:
            address = int(raw_name, 16)
        except ValueError:
            return None
        if not isinstance(self._libs, list) or lib_index >= len(self._libs):
            return None
        return resolve_symbol(self._libs[lib_index], address, self._symbols)


def build_symbol_tables(symbols: dict[str, Any] | None) -> dict[tuple[str, str], SymbolTable]:
    if symbols is None:
        return {}
    string_table = symbols.get("string_table", [])
    tables = {}
    for data in symbols.get("data", []):
        entries = sorted(data.get("symbol_table", []), key=lambda entry: int(entry["rva"]))
        starts = [int(entry["rva"]) for entry in entries]
        key = (
            str(data.get("debug_name", "")),
            str(data.get("code_id", "")).upper(),
        )
        for entry in entries:
            symbol_index = entry.get("symbol")
            if isinstance(symbol_index, int) and symbol_index < len(string_table):
                entry["symbol_name"] = string_table[symbol_index]
        tables[key] = SymbolTable(starts=starts, entries=entries)
    return tables


def resolve_symbol(
    lib: dict[str, Any],
    address: int,
    symbols: dict[tuple[str, str], SymbolTable],
) -> str | None:
    key = (
        str(lib.get("debugName", "")),
        str(lib.get("codeId", "")).upper(),
    )
    table = symbols.get(key)
    if table is None:
        return None
    index = bisect.bisect_right(table.starts, address) - 1
    if index < 0:
        return None
    entry = table.entries[index]
    size = max(1, int(entry.get("size") or 1))
    if address >= int(entry["rva"]) + size:
        return None
    return entry.get("symbol_name")


def format_summary(summary: ProfileSummary) -> str:
    lines = [
        f"thread {summary.thread_name} process={summary.process_name} samples={summary.total_weight}",
        "top inclusive",
    ]
    lines.extend(format_rows(summary.inclusive_rows))
    lines.append("top self")
    lines.extend(format_rows(summary.self_rows))
    return "\n".join(lines)


def format_rows(rows: list[SummaryRow]) -> list[str]:
    return [f"{row.count:8d} {row.percent:6.1f}% {row.name}" for row in rows]


def main() -> None:
    args = parse_args()
    profile = load_profile(args.profile)
    symbols = load_matching_symbols(args.profile) if args.symbols is None else load_symbols(args.symbols)
    summary = summarize_profile(
        profile,
        symbols,
        limit=args.limit,
        pattern=args.pattern,
        thread_name=args.thread,
    )
    print(format_summary(summary))


def load_symbols(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


if __name__ == "__main__":
    main()
