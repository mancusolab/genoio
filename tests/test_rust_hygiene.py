from pathlib import Path

RUST_SOURCE_ROOTS = (
    Path("rust/genoio-core/src"),
    Path("rust/genoio-io/src"),
    Path("rust/genoio-py/src"),
)


def test_rust_production_sources_do_not_unwrap_expect_or_panic():
    offenders: list[str] = []
    for source_root in RUST_SOURCE_ROOTS:
        for path in sorted(source_root.glob("**/*.rs")):
            pending_cfg_test = False
            test_module_depth: int | None = None
            for line_number, line in enumerate(path.read_text().splitlines(), start=1):
                stripped = line.strip()
                if test_module_depth is not None:
                    test_module_depth += line.count("{") - line.count("}")
                    if test_module_depth == 0:
                        test_module_depth = None
                    continue
                if stripped == "#[cfg(test)]":
                    pending_cfg_test = True
                    continue
                if pending_cfg_test and stripped.startswith("mod tests"):
                    test_module_depth = line.count("{") - line.count("}")
                    pending_cfg_test = False
                    continue
                if stripped and not stripped.startswith("#"):
                    pending_cfg_test = False
                if ".unwrap()" in line or ".expect(" in line or "panic!" in line:
                    offenders.append(f"{path}:{line_number}: {stripped}")

    assert offenders == []
