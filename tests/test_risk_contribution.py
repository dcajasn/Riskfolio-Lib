import numpy as np
import pandas as pd
import pytest

import riskfolio as rp

rk = rp.RiskFunctions


@pytest.fixture(scope="module")
def data():
    rng = np.random.default_rng(0)
    returns = pd.DataFrame(
        rng.standard_t(4, (500, 5)) * 0.01 + 0.0003, columns=list("ABCDE")
    )
    w = np.array([0.3, 0.2, 0.2, 0.15, 0.15]).reshape(-1, 1)
    return returns, w, returns.cov().to_numpy()


def _central_differences(fn, w, h):
    n = w.shape[0]
    rc = []
    for i in range(n):
        d = np.zeros((n, 1))
        d[i, 0] = h
        rc.append((fn(w + d) - fn(w - d)) / (2 * h) * w[i, 0])
    return np.array(rc)


@pytest.mark.parametrize("rm", ["EDaR_Rel", "RLDaR_Rel"])
def test_solver_based_relative_drawdown_contributions_use_a_wide_step(data, rm):
    """EDaR_Rel and RLDaR_Rel come from a convex solver, so a 1e-7 step turns their
    finite differences into solver noise. They must use the same 1e-4 step as the
    absolute entropic and relativistic measures, giving contributions that match an
    independent central difference and stay close to the risk itself."""
    returns, w, cov = data
    X = returns.to_numpy()
    if rm == "EDaR_Rel":
        fn = lambda v: rk.EDaR_Rel(X @ v, alpha=0.05)[0]
    else:
        fn = lambda v: rk.RLDaR_Rel(X @ v, alpha=0.05, kappa=0.3)
    rc = np.ravel(rk.Risk_Contribution(w, returns, cov, rm=rm, alpha=0.05, kappa=0.3))
    reference = _central_differences(fn, w, 1e-4)
    np.testing.assert_allclose(rc, reference, rtol=1e-6, atol=1e-8)
    risk = rk.Sharpe_Risk(returns, w, cov, rm=rm, alpha=0.05, kappa=0.3)
    # compounded drawdowns are not homogeneous of degree one, so the Euler sum is
    # only approximately the risk; it must not be off by a large factor
    assert abs(rc.sum() - risk) / risk < 0.1


@pytest.mark.parametrize("rm", ["MV", "MAD", "CVaR", "MDD", "CDaR", "UCI"])
def test_euler_identity_for_homogeneous_measures(data, rm):
    returns, w, cov = data
    rc = rk.Risk_Contribution(w, returns, cov, rm=rm, alpha=0.05)
    risk = rk.Sharpe_Risk(returns, w, cov, rm=rm, alpha=0.05)
    assert np.sum(rc) == pytest.approx(risk, rel=1e-5)


def test_unknown_risk_measure_raises_value_error(data):
    returns, w, cov = data
    with pytest.raises(ValueError, match="valid risk measure"):
        rk.Sharpe_Risk(returns, w, cov, rm="NOT_A_MEASURE")
    with pytest.raises(ValueError, match="valid risk measure"):
        rk.Risk_Contribution(w, returns, cov, rm="NOT_A_MEASURE")
