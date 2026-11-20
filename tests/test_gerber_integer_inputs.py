import numpy as np
import pandas as pd
import pytest

from riskfolio.src.GerberStatistic import gerber_cov_stat2


@pytest.mark.parametrize(
    "dtype,n", [(np.int8, 200), (np.int16, 40000), (np.int32, 200)]
)
@pytest.mark.parametrize("dataframe", [False, True])
def test_gerber_integer_counts_do_not_overflow(dtype, n, dataframe):
    returns = np.tile([[-1, 1], [1, -1]], (n // 2, 1)).astype(dtype)
    expected = np.array([[1.0, -1.0], [-1.0, 1.0]])
    if dataframe:
        returns = pd.DataFrame(returns, columns=["asset_a", "asset_b"])

    with np.errstate(invalid="raise", divide="raise"):
        result = gerber_cov_stat2(returns)

    np.testing.assert_allclose(result, expected)
    if dataframe:
        assert result.index.tolist() == ["asset_a", "asset_b"]
        assert result.columns.tolist() == ["asset_a", "asset_b"]


def test_gerber_fractional_returns_keep_threshold_precision():
    returns = np.array([[0.01, -0.02], [-0.01, 0.02]] * 100)
    expected = np.array([[0.0001, -0.0002], [-0.0002, 0.0004]])
    np.testing.assert_allclose(gerber_cov_stat2(returns), expected)
