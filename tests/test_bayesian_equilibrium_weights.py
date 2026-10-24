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


@pytest.mark.parametrize("rf", [0.0, 0.004])
@pytest.mark.parametrize("alpha", [0.0, 0.002, -0.001])
def test_bayesian_one_asset_matches_scalar_posterior(rf, alpha):
    # Scalar Gaussian conditioning gives an independent reference without
    # repeating the implementation's nested matrix inversions.
    factors = pd.DataFrame({"F": [-0.02, 0.01, 0.04, 0.03, -0.01]})
    residuals = np.array([0.003, -0.002, 0.001, -0.004, 0.002])
    loading, delta, view = 1.2, 2.0, 0.01
    assets = pd.DataFrame(
        {"A": alpha + loading * factors["F"].to_numpy() + residuals}
    )
    factor_variance = factors["F"].var(ddof=1)
    tau = 1 / len(factors)
    posterior_factor_mean = (tau * (factors["F"].mean() - rf) + view) / (1 + tau)
    posterior_factor_variance = tau * factor_variance / (1 + tau)
    expected_variance = (
        loading**2 * (factor_variance + posterior_factor_variance)
        + residuals.var(ddof=1)
    )
    expected_excess_return = alpha + loading * posterior_factor_mean

    mu, cov, weights = black_litterman_bayesian(
        assets,
        factors,
        pd.DataFrame([[alpha, loading]]),
        np.array([[1.0]]),
        np.array([[view]]),
        delta=delta,
        rf=rf,
        const=True,
    )
    np.testing.assert_allclose(mu.iloc[0, 0], expected_excess_return + rf)
    np.testing.assert_allclose(cov.iloc[0, 0], expected_variance)
    # Maximizer of w * excess_return - delta/2 * variance * w**2.
    np.testing.assert_allclose(
        weights.iloc[0, 0], expected_excess_return / (delta * expected_variance)
    )
