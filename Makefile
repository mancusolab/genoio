.PHONY: build-dev build-release check docs pyright rust-check rust-doc rust-fmt rust-test test verify

PYTHON ?= python
CARGO ?= cargo
MKDOCS ?= mkdocs
PYTEST ?= pytest
MATURIN ?= $(PYTHON) -m maturin
# macOS builds need explicit tool names in some managed Python environments.
RUST_ENV ?= env CC=clang AR=ar

build-dev:
	$(RUST_ENV) $(MATURIN) develop

build-release:
	$(RUST_ENV) $(MATURIN) develop --release

test:
	$(PYTEST) -q

pyright:
	pyright src tests scripts

docs:
	$(MKDOCS) build --strict

rust-fmt:
	$(RUST_ENV) $(CARGO) fmt --manifest-path rust/Cargo.toml --all -- --check

rust-check:
	$(RUST_ENV) $(CARGO) clippy --manifest-path rust/Cargo.toml --workspace --all-targets -- -D warnings

rust-test:
	$(RUST_ENV) $(CARGO) test --manifest-path rust/Cargo.toml --workspace

rust-doc:
	$(RUST_ENV) RUSTDOCFLAGS=-Dwarnings $(CARGO) doc --manifest-path rust/Cargo.toml --workspace --no-deps

check: rust-fmt rust-check rust-test pyright test docs

verify: build-dev check rust-doc
