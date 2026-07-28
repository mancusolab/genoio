.DEFAULT_GOAL := help
SHELL := /bin/bash

.PHONY: build build-dev build-release build-wheel check clean docs fresh help lock requirements ruff-check ruff-fmt rust-audit rust-check rust-doc rust-fmt rust-test test ty venv verify

VENV ?= .venv
SYSTEM_PYTHON ?= python3
UV ?= uv

ifeq ($(OS),Windows_NT)
	VENV_BIN := $(VENV)/Scripts
else
	VENV_BIN := $(VENV)/bin
endif

PYTHON ?= $(VENV_BIN)/python
CARGO ?= cargo
ZENSICAL ?= $(VENV_BIN)/zensical
TY ?= $(VENV_BIN)/ty
PYTEST ?= $(VENV_BIN)/pytest
RUFF ?= $(VENV_BIN)/ruff
MATURIN ?= env -u CONDA_PREFIX VIRTUAL_ENV=$(abspath $(VENV)) $(PYTHON) -m maturin
DIST_DIR ?= dist
CODESIGN_ALLOCATE ?= /usr/bin/codesign_allocate
ARGS ?=
# macOS builds need explicit tool names in some managed Python environments.
RUST_ENV ?= env CC=clang AR=ar
REPAIR_ENV ?= env CODESIGN_ALLOCATE=$(CODESIGN_ALLOCATE) CC=clang AR=ar

venv: requirements  ## Create and sync the local Python virtual environment

lock:  ## Update uv.lock from pyproject.toml
	$(UV) lock

requirements: uv.lock  ## Sync development and documentation dependencies
	$(RUST_ENV) UV_PROJECT_ENVIRONMENT=$(VENV) $(UV) sync --frozen --all-extras --no-install-project --python $(SYSTEM_PYTHON)

build: build-dev  ## Compile and install the extension for development

build-dev: requirements  ## Compile and install the extension for development
	$(RUST_ENV) $(MATURIN) develop $(ARGS)

build-release: requirements  ## Compile and install an optimized development extension
	$(RUST_ENV) $(MATURIN) develop --release $(ARGS)

build-wheel: requirements  ## Build a repaired redistributable wheel
	$(REPAIR_ENV) $(MATURIN) build --release --auditwheel=repair -o $(DIST_DIR) $(ARGS)
	$(RUST_ENV) $(MATURIN) develop

test: build-dev  ## Run Python tests
	$(PYTEST) -p no:capture -q

ty: build-dev  ## Run Python type checks
	$(TY) check src tests scripts

ruff-check: requirements  ## Run Python lint checks
	$(RUFF) check src tests scripts

ruff-fmt: requirements  ## Check Python formatting
	$(RUFF) format --check src tests scripts

docs: build-dev  ## Build documentation with strict checks
	$(ZENSICAL) build --strict
	rm -rf site

rust-fmt:  ## Check Rust formatting
	$(RUST_ENV) $(CARGO) fmt --manifest-path rust/Cargo.toml --all -- --check

rust-check:  ## Run Rust clippy checks
	$(RUST_ENV) $(CARGO) clippy --manifest-path rust/Cargo.toml --workspace --all-targets -- -D warnings

rust-audit:  ## Audit Rust dependencies with cargo-audit
	cd rust && $(CARGO) audit

rust-test:  ## Run Rust tests
	$(RUST_ENV) $(CARGO) test --manifest-path rust/Cargo.toml --workspace

rust-doc:  ## Build Rust docs with warnings as errors
	$(RUST_ENV) RUSTDOCFLAGS=-Dwarnings $(CARGO) doc --manifest-path rust/Cargo.toml --workspace --no-deps

check: rust-fmt rust-check rust-test ruff-check ruff-fmt ty test docs  ## Run the standard validation suite

verify: build-dev check rust-doc  ## Build the extension and run all validation checks

fresh: clean build-dev  ## Recreate the local environment and rebuild the extension

clean:  ## Remove local build and test artifacts
	rm -rf build dist dist-repaired dist-repaired-ok site .mypy_cache .pytest_cache .ruff_cache
	find . -type d -name __pycache__ -prune -exec rm -rf {} +

help:  ## Display this help screen
	@echo "Available commands:"
	@grep -E '^[a-zA-Z0-9_.-]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-18s\033[0m %s\n", $$1, $$2}' | sort
	@echo
	@echo 'Run make lock after changing Python dependencies.'
	@echo 'Build targets accept ARGS, for example: make build-wheel ARGS="--sdist"'
