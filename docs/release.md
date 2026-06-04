# Release

This page records the local release checks for maintainers. It is intentionally
short: the public documentation covers usage, while this page keeps packaging
and documentation deployment steps in one place.

## Verify the tree

Run the full local verification suite from a clean worktree:

```bash
make verify
```

This builds the Rust extension, checks Rust formatting and lints, runs Rust and
Python tests, runs Pyright, builds the MkDocs site with strict checks, and builds
Rust docs with warnings treated as errors.

## Build a wheel

Build a repaired wheel before testing installation:

```bash
make build-wheel
```

On macOS this uses `--auditwheel=repair` and an explicit `CODESIGN_ALLOCATE`
path so external dynamic libraries are bundled into the wheel instead of
pointing at a local environment.

## Smoke-test the artifact

Install the repaired wheel into a fresh environment and test each supported
format with a small fixture:

```bash
python -m venv /tmp/genoio-wheel-test
/tmp/genoio-wheel-test/bin/python -m pip install dist/genoio-*.whl
/tmp/genoio-wheel-test/bin/python -c "import genoio; print(genoio.__version__)"
```

For beta releases, smoke-test VCF, PLINK1, and PLINK2 reads before publishing.

## Deploy documentation

After GitHub Pages is configured for the repository, deploy the docs from the
release commit:

```bash
zensical build --strict
```

Publish the generated `site/` directory with the repository's GitHub Pages
workflow or another static-site deployment step.

The configured documentation URL is
`https://mancusolab.github.io/genoio`.
