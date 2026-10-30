import numpy as np
import pandas as pd

import riskfolio as rp


def test_cov_returns_keeps_dataframe_columns():
    cov = pd.DataFrame(
        [[0.04, 0.01], [0.01, 0.09]], index=["A", "B"], columns=["A", "B"]
    )
    returns = rp.cov_returns(cov)
    assert isinstance(returns, pd.DataFrame)
    assert returns.columns.tolist() == ["A", "B"]
    np.testing.assert_allclose(np.cov(returns.to_numpy().T), cov.to_numpy())


def test_cov_returns_ndarray_input_returns_ndarray():
    cov = np.array([[0.04, 0.01], [0.01, 0.09]])
    returns = rp.cov_returns(cov)
    assert isinstance(returns, np.ndarray)
    np.testing.assert_allclose(np.cov(returns.T), cov)
