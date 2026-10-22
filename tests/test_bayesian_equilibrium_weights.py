import numpy as np
import pandas as pd
import pytest

from riskfolio.src.ParamsEstimation import black_litterman_bayesian


@pytest.mark.parametrize("rf", [0.0, 0.004])
@pytest.mark.parametrize("const", [False, True])
def test_bayesian_weights_satisfy_excess_return_first_order_condition(rf, const):
    rng = np.random.default_rng(13)
    factors = pd.DataFrame(rng.normal(size=(100, 2)) * 0.03 + 0.002)
    loadings = np.array([[1.0, 0.2], [0.8, -0.3], [1.2, 0.5]])
    assets = pd.DataFrame(
        factors.to_numpy() @ loadings.T + rng.normal(size=(100, 3)) * 0.02,
        columns=list("ABC"),
    )
    if const:
        loadings = np.column_stack(([0.001, 0.002, -0.001], loadings))
    mu, cov, weights = black_litterman_bayesian(
        assets,
        factors,
        pd.DataFrame(loadings),
        np.array([[1.0, 0.0]]),
        np.array([[0.01]]),
        delta=2,
        rf=rf,
        const=const,
    )
    # Stationarity of w'(mu-rf) - delta/2 w'cov w, including any estimated alpha.
    np.testing.assert_allclose(
        2 * cov.to_numpy() @ weights.to_numpy(),
        mu.to_numpy().T - rf,
        rtol=1e-10,
        atol=1e-12,
    )
