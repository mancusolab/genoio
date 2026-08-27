# Release

Publishing a GitHub Release starts two workflows from the release tag:

- **Release** builds and smoke-tests eight wheels and one source distribution,
  attaches them to the GitHub Release, and publishes them to PyPI with trusted
  publishing.
- **Documentation** builds the versioned package API and deploys the site to
  GitHub Pages.

Both workflows reject a tag that does not match the package version. A manual
workflow run requires an explicit release tag and performs the same publishing
or deployment operation; it is not a dry run.

## Prepare the release commit

Start from a clean `main` branch with passing CI and Build workflows. Update the
version in these files:

- `pyproject.toml`
- `rust/Cargo.toml`
- `src/genoio/__init__.py` (the source-tree fallback)

Refresh both lockfiles and confirm that only the genoio package versions changed:

```bash
uv lock
cargo check --manifest-path rust/Cargo.toml --workspace
git diff -- uv.lock rust/Cargo.lock
```

Run the complete verification suite:

```bash
make verify
```

This syncs the locked Python environment, builds the Rust extension, checks Rust
formatting and lints, runs Rust and Python tests, runs Ty, builds the Zensical
site with strict checks, and builds Rust docs with warnings treated as errors.

Commit and push the version update. Wait for both CI and Build to pass on that
commit before creating the release tag.

## Verify a wheel locally

Build a repaired wheel:

```bash
make build-wheel
```

On macOS, this uses `--auditwheel=repair` and an explicit `CODESIGN_ALLOCATE`
path so external dynamic libraries are bundled instead of pointing at the local
environment. GitHub Actions builds Linux wheels against `manylinux_2_28` and
checks PyPI compatibility for both Linux and macOS wheels.

Install the repaired wheel into a fresh environment and run the smoke test:

```bash
uv venv /tmp/genoio-wheel-test
uv pip install --python /tmp/genoio-wheel-test/bin/python dist/genoio-*.whl
/tmp/genoio-wheel-test/bin/python scripts/wheel_smoke.py
```

## Publish

Create and push a tag that points to the verified release commit, then publish
the GitHub Release:

```bash
git tag v0.4.1
git push origin v0.4.1
gh release create v0.4.1 --verify-tag --generate-notes
```

Replace `v0.4.1` with the intended version. Publishing the GitHub Release starts
the Release and Documentation workflows. The Release workflow publishes to the
`pypi` GitHub environment, which must remain configured as a trusted publisher
for the `genoio` PyPI project.

Monitor both workflows and verify the uploaded version:

```bash
gh run list --workflow release.yml --limit 1
gh run list --workflow deploydocs.yml --limit 1
python -m pip index versions genoio
```

The documentation site is `https://mancusolab.github.io/genoio`.
