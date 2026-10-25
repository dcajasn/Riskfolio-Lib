import numpy as np
import pandas as pd
import pytest

import riskfolio as rp
import riskfolio.src.ParamsEstimation as pe

RNG = np.random.default_rng(0)
X = pd.DataFrame(RNG.normal(0.0005, 0.01, (500, 5)), columns=list("ABCDE"))


@pytest.mark.parametrize("method", ["JS", "BS", "BOP"])
@pytest.mark.parametrize("target", ["b1", "b2", "b3"])
def test_shrinkage_mean_is_real(method, target):
    # np.linalg.eigvals returns complex128 for symmetric matrices on recent
    # numpy, which used to turn the James-Stein mean into a complex vector.
    mu = pe.mean_vector(X, method=method, target=target)
    assert list(mu.columns) == list(X.columns)
    assert mu.shape == (1, X.shape[1])
    assert not np.iscomplexobj(mu.to_numpy())
    assert np.isfinite(mu.to_numpy()).all()


def test_james_stein_matches_closed_form():
    T, n = X.shape
    m = X.mean().to_numpy().reshape(-1, 1)
    ev = np.linalg.eigvalsh(np.cov(X, rowvar=False))
    b = np.full((n, 1), m.mean())
    a = (n * ev.mean() - 2 * ev.max()) / T / ((m - b).T @ (m - b)).item()
    expected = (1 - a) * m + a * b

    mu = pe.mean_vector(X, method="JS", target="b1")
    np.testing.assert_allclose(mu.to_numpy().ravel(), expected.ravel(), rtol=1e-12)


def test_js_mean_can_be_optimized():
    # A complex mu made cvxpy fail with "Inequality constraints cannot be complex".
    port = rp.Portfolio(returns=X)
    port.assets_stats(method_mu="JS", method_cov="hist")
    w = port.optimization(model="Classic", rm="MV", obj="Sharpe")
    assert w is not None
    assert w["weights"].sum() == pytest.approx(1.0)
