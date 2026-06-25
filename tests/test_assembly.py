# pattern: Imperative Shell

import numpy as np
import pytest

from genoio._assembly import dense_array_from_rust, sparse_matrix_from_rust


def test_dense_array_from_rust_assembles_variant_major_payload_as_strided_view():
    array = dense_array_from_rust(
        values=np.array([0.0, 1.0, 2.0, 3.0], dtype=np.float32),
        shape=(2, 2),
        values_layout="variant_major",
        dtype=np.dtype("float32"),
    )

    np.testing.assert_array_equal(
        array,
        np.array([[0.0, 2.0], [1.0, 3.0]], dtype=np.float32),
    )
    assert array.flags.writeable
    assert not array.flags.c_contiguous


def test_dense_array_from_rust_rejects_unknown_dense_layout():
    with pytest.raises(AssertionError, match="dense value layout"):
        dense_array_from_rust(
            values=np.array([0.0], dtype=np.float32),
            shape=(1, 1),
            dtype=np.dtype("float32"),
            values_layout="columnar",
        )


def test_dense_array_from_rust_rejects_values_that_do_not_match_shape():
    with pytest.raises(AssertionError, match="does not match shape"):
        dense_array_from_rust(
            values=np.array([0.0, 1.0], dtype=np.float32),
            shape=(1, 1),
            dtype=np.dtype("float32"),
        )


def test_sparse_matrix_from_rust_preserves_int32_indices():
    matrix = sparse_matrix_from_rust(
        data=np.array([1.0, 2.0], dtype=np.float32),
        indices=np.array([1, 2], dtype=np.int32),
        indptr=np.array([0, 2], dtype=np.int32),
        shape=(3, 1),
        dtype=np.dtype("float32"),
        sparse_format="csc",
    )

    assert matrix.indices.dtype == np.dtype("int32")
    assert matrix.indptr.dtype == np.dtype("int32")


def test_sparse_matrix_from_rust_rejects_unsupported_index_dtype():
    with pytest.raises(AssertionError, match="sparse indices dtype"):
        sparse_matrix_from_rust(
            data=np.array([1.0], dtype=np.float32),
            indices=np.array([0], dtype=np.int64),
            indptr=np.array([0, 1], dtype=np.int32),
            shape=(1, 1),
            dtype=np.dtype("float32"),
            sparse_format="csc",
        )
