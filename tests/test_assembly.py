# pattern: Imperative Shell

import numpy as np
import pytest

import genoio._assembly as assembly
from genoio._assembly import _impute_missing_by_variant, dense_array_from_rust


def test_impute_missing_by_variant_returns_original_array_when_no_values_are_missing():
    array = np.array([[0.0, 1.0], [2.0, 0.0]], dtype=np.float32)
    mask = np.zeros(array.shape, dtype=bool)

    imputed = _impute_missing_by_variant(array, mask)

    assert imputed is array


def test_dense_array_from_rust_assembles_variant_major_payload_as_strided_view():
    array = dense_array_from_rust(
        values=np.array([0.0, 1.0, 2.0, 3.0], dtype=np.float32),
        shape=(2, 2),
        missing_mask=np.array([False, True, False, False], dtype=bool),
        values_layout="variant_major",
        missing="nan",
        dtype=np.dtype("float32"),
    )

    np.testing.assert_array_equal(
        array,
        np.array([[0.0, 2.0], [np.nan, 3.0]], dtype=np.float32),
    )
    assert array.flags.writeable
    assert not array.flags.c_contiguous


def test_dense_array_from_rust_imputes_from_variant_major_missing_indices():
    array = dense_array_from_rust(
        values=np.array([0.0, 1.0, 2.0, 3.0], dtype=np.float32),
        shape=(2, 2),
        missing="impute",
        dtype=np.dtype("float32"),
        values_layout="variant_major",
        missing_indices=np.array([1], dtype=np.int64),
    )

    np.testing.assert_array_equal(
        array,
        np.array([[0.0, 2.0], [0.0, 3.0]], dtype=np.float32),
    )


def test_dense_array_from_rust_imputes_sparse_indices_without_dense_mask(monkeypatch):
    def fail_dense_mask(*args, **kwargs):
        raise AssertionError("sparse index imputation should not build a dense mask")

    monkeypatch.setattr(assembly, "_dense_missing_mask_from_flat", fail_dense_mask)

    array = dense_array_from_rust(
        values=np.array([2.0, 10.0, 100.0, 20.0, 6.0, 30.0], dtype=np.float32),
        shape=(3, 2),
        missing="impute",
        dtype=np.dtype("float32"),
        values_layout="sample_major",
        missing_indices=np.array([2], dtype=np.int64),
    )

    np.testing.assert_array_equal(
        array,
        np.array([[2.0, 10.0], [4.0, 20.0], [6.0, 30.0]], dtype=np.float32),
    )


def test_dense_array_from_rust_imputes_sparse_indices_in_place():
    values = np.array([2.0, 10.0, 100.0, 20.0, 6.0, 30.0], dtype=np.float32)

    array = dense_array_from_rust(
        values=values,
        shape=(3, 2),
        missing="impute",
        dtype=np.dtype("float32"),
        values_layout="sample_major",
        missing_indices=np.array([2], dtype=np.int64),
    )

    assert np.shares_memory(array, values)
    np.testing.assert_array_equal(
        values,
        np.array([2.0, 10.0, 4.0, 20.0, 6.0, 30.0], dtype=np.float32),
    )


def test_dense_array_from_rust_accepts_empty_missing_index_sequence():
    values = np.array([1.0], dtype=np.float32)

    array = dense_array_from_rust(
        values=values,
        shape=(1, 1),
        missing="impute",
        dtype=np.dtype("float32"),
        missing_indices=[],
    )

    assert np.shares_memory(array, values)
    np.testing.assert_array_equal(array, np.array([[1.0]], dtype=np.float32))


def test_dense_array_from_rust_rejects_all_missing_variant_from_sparse_indices():
    with pytest.raises(assembly.MissingDataError, match="all-missing variant"):
        dense_array_from_rust(
            values=np.array([100.0, 10.0, 200.0, 20.0], dtype=np.float32),
            shape=(2, 2),
            missing="impute",
            dtype=np.dtype("float32"),
            values_layout="sample_major",
            missing_indices=np.array([0, 2], dtype=np.int64),
        )


def test_dense_array_from_rust_rejects_non_integer_missing_indices():
    with pytest.raises(AssertionError, match="integer"):
        dense_array_from_rust(
            values=np.array([0.0], dtype=np.float32),
            shape=(1, 1),
            missing="impute",
            dtype=np.dtype("float32"),
            missing_indices=np.array([0.0], dtype=np.float32),
        )


def test_dense_array_from_rust_rejects_conflicting_missing_payloads():
    with pytest.raises(AssertionError, match="mask or indices"):
        dense_array_from_rust(
            values=np.array([0.0], dtype=np.float32),
            shape=(1, 1),
            missing="impute",
            dtype=np.dtype("float32"),
            missing_mask=np.array([True], dtype=bool),
            missing_indices=np.array([0], dtype=np.int64),
        )
