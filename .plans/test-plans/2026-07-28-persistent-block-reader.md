# Persistent block reader human test plan

**Reviewed implementation:** `0b64e5be1109447c1563ac5756305c33ca565737`

**Status:** Supplemental environmental validation. Automated tests cover every
acceptance criterion; these scenarios exercise production-scale files and
operating-system behavior that small deterministic fixtures don't represent.

## Automated coverage summary

- All 20 leaf acceptance criteria have passing automated evidence.
- All 27 normalized requirement IDs are present in test names.
- The final verification run passed 442 Rust tests and 450 Python tests.
- Two Python tests skipped because optional local metadata fixtures were absent.
  Neither skip covers persistent block reading.

## Preconditions

- Run commands from the repository worktree.
- Build the current native extension with `make build-dev`.
- Use real, nonempty files supplied through absolute paths.
- Install `lsof` for the handle-release scenario.
- Supply `.tbi` or `.csi` beside an indexed VCF and `.bgi` beside an indexed
  BGEN.

Record the input file versions, command output, and observed result for each
scenario.

## 1. Observe memory use on a large VCF

Run the scan twice, once with `PBR_BLOCK_SIZE=1024` and once with
`PBR_BLOCK_SIZE=8192`.

```bash
PBR_LARGE_VCF=/absolute/path/cohort.vcf.gz \
PBR_BLOCK_SIZE=1024 \
/usr/bin/time -l .venv/bin/python - <<'PY'
import os

import genoio

size = int(os.environ["PBR_BLOCK_SIZE"])
dataset = genoio.vcf(os.environ["PBR_LARGE_VCF"])
variants = 0
blocks = 0

for matrix, metadata in dataset.iter_blocks(
    size=size,
    return_variants=True,
):
    assert 0 < matrix.shape[1] <= size
    assert matrix.shape[1] == metadata.height
    variants += matrix.shape[1]
    blocks += 1

print({"blocks": blocks, "variants": variants})
PY
```

Expected results:

- The scan completes without constructing a whole-dataset genotype matrix.
- Every block contains at most the requested number of variants.
- Each metadata frame has one row per matrix column.
- Maximum RSS changes mainly with block size and decoder scratch, not total
  variant count. Timing is descriptive and isn't a pass/fail threshold.

Related criteria: AC1.2, AC1.3, AC4.1, AC4.2.

## 2. Inspect handle release after early close

Start the reader in terminal 1:

```bash
PBR_LARGE_VCF=/absolute/path/cohort.vcf.gz .venv/bin/python - <<'PY'
import os

import genoio

iterator = genoio.vcf(os.environ["PBR_LARGE_VCF"]).iter_blocks(1024)
next(iterator)
print("PID", os.getpid(), flush=True)
input("Inspect open handles, then press Enter to close: ")
iterator.close()
input("Inspect handles again, then press Enter to exit: ")
PY
```

In terminal 2, replace `PID` with the printed process ID:

```bash
lsof -p PID | rg 'cohort\.vcf\.gz|\.tbi|\.csi'
```

Expected results:

- Before `close()`, the source and resolved index each appear at most once.
- After `close()`, those handles are absent.
- The iterator releases its handles without being exhausted.

Related criterion: AC3.4.

## 3. Compare real indexed VCF and BGEN reads

Choose a region containing at least one retained variant in both files.

```bash
PBR_INDEXED_VCF=/absolute/path/cohort.vcf.gz \
PBR_INDEXED_BGEN=/absolute/path/cohort.bgen \
PBR_REGION=1:100000-200000 \
.venv/bin/python - <<'PY'
import os

import numpy as np

import genoio

region = genoio.region(os.environ["PBR_REGION"])
cases = [
    ("VCF", genoio.vcf(os.environ["PBR_INDEXED_VCF"]), {}),
    ("BGEN", genoio.bgen(os.environ["PBR_INDEXED_BGEN"]), {"dosage": "dosage"}),
]

for label, dataset, options in cases:
    full, full_variants = dataset.read(
        variants=region,
        return_variants=True,
        **options,
    )
    blocks = list(
        dataset.iter_blocks(
            4096,
            variants=region,
            return_variants=True,
            **options,
        )
    )
    assert blocks, f"choose a nonempty region for {label}"

    combined = np.concatenate([matrix for matrix, _ in blocks], axis=1)
    rows = [row for _, variants in blocks for row in variants.rows()]
    np.testing.assert_allclose(combined, full, equal_nan=True)
    assert rows == full_variants.rows()
    print(label, combined.shape)
PY
```

Expected results:

- Concatenated blocks equal the corresponding whole read.
- Variant metadata remains in source order.
- Index or chunk boundaries don't duplicate or omit variants.

Related criteria: AC1.2, AC2.3.

## 4. Exercise production PLINK2 storage modes

Run the command against representative fixed-width, variable-width, phased,
dosage, and LD-compressed PGEN files. Set `PBR_READ_OPTIONS` for the mode under
test, for example:

- `{}`
- `{"sparse": "csc"}`
- `{"dosage": "dosage"}`
- `{"kind": "haplo", "dosage": "dosage"}`

```bash
PBR_PFILE=/absolute/path/prefix \
PBR_READ_OPTIONS='{"dosage": "dosage"}' \
.venv/bin/python - <<'PY'
import json
import os

import numpy as np
from scipy import sparse

import genoio

options = json.loads(os.environ["PBR_READ_OPTIONS"])
dataset = genoio.pfile(os.environ["PBR_PFILE"])
full = dataset.read(**options)
blocks = list(dataset.iter_blocks(1024, **options))
assert blocks

if sparse.issparse(full):
    combined = sparse.hstack(
        blocks,
        format=options.get("sparse", "csc"),
    )
    np.testing.assert_allclose(
        combined.toarray(),
        full.toarray(),
        equal_nan=True,
    )
else:
    combined = np.concatenate(blocks, axis=1)
    np.testing.assert_allclose(combined, full, equal_nan=True)

print(combined.shape)
PY
```

Expected results:

- Each supported production mode matches `Dataset.read()`.
- Dependencies that cross block boundaries don't change decoded values.
- Unsupported modes retain their documented public exception class.

Related criteria: AC2.1, AC2.4, AC2.7.

## Result record

| Scenario | Input/corpus | Result | Notes |
|---|---|---|---|
| Large-file memory |  |  |  |
| Early-close handles |  |  |  |
| Indexed VCF/BGEN |  |  |  |
| PLINK2 storage modes |  |  |  |
