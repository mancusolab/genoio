# pattern: Imperative Shell

import numpy as np

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
