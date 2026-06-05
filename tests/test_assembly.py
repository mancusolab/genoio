import numpy as np

from genoio._assembly import _impute_missing_by_variant


def test_impute_missing_by_variant_returns_original_array_when_no_values_are_missing():
    array = np.array([[0.0, 1.0], [2.0, 0.0]], dtype=np.float32)
    mask = np.zeros(array.shape, dtype=bool)

    imputed = _impute_missing_by_variant(array, mask)

    assert imputed is array
