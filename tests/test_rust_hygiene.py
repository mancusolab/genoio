import re
from pathlib import Path


def test_internal_contract_errors_map_to_private_internal_error():
    source = Path("rust/genoio-py/src/lib.rs").read_text()
    branch_match = re.search(r"GenoioError::InternalContract\s*\{ \.\. \}\s*=>.*", source)

    assert branch_match is not None
    branch = branch_match.group(0)
    assert "RustInternalError::new_err" in branch
    assert "PyRuntimeError" not in branch
