.DEFAULT_GOAL := help
SHELL := /bin/bash

.PHONY: build build-dev build-release build-wheel check clean docs fresh help pyright requirements rust-check rust-doc rust-fmt rust-test test venv verify

VENV ?= .venv
SYSTEM_PYTHON ?= python3

ifeq ($(OS),Windows_NT)
	VENV_BIN := $(VENV)/Scripts
else
	VENV_BIN := $(VENV)/bin
endif

PYTHON ?= $(VENV_BIN)/python
CARGO ?= cargo
MKDOCS ?= $(VENV_BIN)/mkdocs
PYRIGHT ?= $(VENV_BIN)/pyright
PYTEST ?= $(VENV_BIN)/pytest
MATURIN ?= $(PYTHON) -m maturin
DIST_DIR ?= dist
CODESIGN_ALLOCATE ?= /usr/bin/codesign_allocate
ARGS ?=
# macOS builds need explicit tool names in some managed Python environments.
RUST_ENV ?= env CC=clang AR=ar
REPAIR_ENV ?= env CODESIGN_ALLOCATE=$(CODESIGN_ALLOCATE) CC=clang AR=ar

venv:  ## Create the local Python virtual environment
	$(SYSTEM_PYTHON) -m venv $(VENV)

requirements: venv  ## Install development and documentation dependencies
	$(PYTHON) -m pip install --upgrade pip
	$(RUST_ENV) $(PYTHON) -m pip install -e ".[dev,docs]"

build: build-dev  ## Compile and install the extension for development

build-dev: requirements  ## Compile and install the extension for development
	$(RUST_ENV) $(MATURIN) develop $(ARGS)

build-release: requirements  ## Compile and install an optimized development extension
	$(RUST_ENV) $(MATURIN) develop --release $(ARGS)

build-wheel: requirements  ## Build a repaired redistributable wheel
	$(REPAIR_ENV) $(MATURIN) build --release --auditwheel=repair -o $(DIST_DIR) $(ARGS)

test: requirements  ## Run Python tests
	$(PYTEST) -q

pyright: requirements  ## Run Python type checks
	$(PYRIGHT) src tests scripts

docs: requirements  ## Build documentation with strict checks
	$(MKDOCS) build --strict
	rm -rf site

rust-fmt:  ## Check Rust formatting
	$(RUST_ENV) $(CARGO) fmt --manifest-path rust/Cargo.toml --all -- --check

rust-check:  ## Run Rust clippy checks
	$(RUST_ENV) $(CARGO) clippy --manifest-path rust/Cargo.toml --workspace --all-targets -- -D warnings

rust-test:  ## Run Rust tests
	$(RUST_ENV) $(CARGO) test --manifest-path rust/Cargo.toml --workspace

rust-doc:  ## Build Rust docs with warnings as errors
	$(RUST_ENV) RUSTDOCFLAGS=-Dwarnings $(CARGO) doc --manifest-path rust/Cargo.toml --workspace --no-deps

check: rust-fmt rust-check rust-test pyright test docs  ## Run the standard validation suite

verify: build-dev check rust-doc  ## Build the extension and run all validation checks

fresh: clean requirements build-dev  ## Recreate the local environment and rebuild the extension

clean:  ## Remove local build and test artifacts
	rm -rf build dist dist-repaired dist-repaired-ok site .mypy_cache .pytest_cache .ruff_cache
	find . -type d -name __pycache__ -prune -exec rm -rf {} +

help:  ## Display this help screen
	@echo "Available commands:"
	@grep -E '^[a-zA-Z0-9_.-]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-18s\033[0m %s\n", $$1, $$2}' | sort
	@echo
	@echo 'Build targets accept ARGS, for example: make build-wheel ARGS="--sdist"'
