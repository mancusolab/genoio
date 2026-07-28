# pattern: Imperative Shell

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest
from script_loader import load_benchmark_script

wheel_smoke = load_benchmark_script("wheel_smoke")


def test_make_test_target_disables_pytest_capture_plugin() -> None:
    result = subprocess.run(
        ["make", "-n", "test", "PYTEST=pytest"],
        check=True,
        capture_output=True,
        text=True,
    )
    pytest_commands = [line for line in result.stdout.splitlines() if line.startswith("pytest ")]

    assert pytest_commands == ["pytest -p no:capture -q"]


def test_wheel_smoke_rejects_install_without_typing_marker(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    monkeypatch.setattr(wheel_smoke, "files", lambda package: tmp_path, raising=False)

    with pytest.raises(AssertionError, match="py\\.typed"):
        wheel_smoke.main()
