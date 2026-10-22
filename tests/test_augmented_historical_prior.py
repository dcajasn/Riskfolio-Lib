import numpy as np
import pandas as pd
import pytest

from riskfolio.src.ParamsEstimation import augmented_black_litterman


@pytest.mark.parametrize("views", ["assets", "factors", "both"])
@pytest.mark.parametrize("eq", [False, True])
def test_augmented_prior_matches_gaussian_conditioning(views, eq):
    rng = np.random.default_rng(13)
    factors = pd.DataFrame(rng.normal(size=(100, 2)) * 0.03 + 0.002)
    b = np.array([[1.0, 0.2], [0.8, -0.3], [1.2, 0.5]])
    assets = pd.DataFrame(
        factors.to_numpy() @ b.T + rng.normal(size=(100, 3)) * 0.02,
        columns=list("ABC"),
    )
    w = pd.DataFrame([0.4, 0.35, 0.25], index=assets.columns)
    rf, delta = 0.004, 2.0
    pa, pf = np.array([[1.0, -1.0, 0.0]]), np.array([[1.0, 0.0]])
    qa, qf = np.array([[0.002]]), np.array([[0.01]])
    kwargs = {}
    rows, targets = [], []
    if views in ["assets", "both"]:
        kwargs.update(P=pa, Q=qa)
        rows.append(np.hstack((pa, np.zeros((1, 2)))))
        targets.append(qa)
    if views in ["factors", "both"]:
        kwargs.update(P_f=pf, Q_f=qf)
        rows.append(np.hstack((np.zeros((1, 3)), pf)))
        targets.append(qf)
    mu, cov, weights = augmented_black_litterman(
        assets,
        w,
        factors,
        pd.DataFrame(b),
        delta=delta,
        rf=rf,
        eq=eq,
        const=False,
        **kwargs,
    )
    s, sf = assets.cov().to_numpy(), factors.cov().to_numpy()
    sa = np.block([[s, b @ sf], [sf @ b.T, sf]])
    prior = (
        delta * np.vstack((s, sf @ b.T)) @ w.to_numpy()
        if eq
        else np.concatenate((assets.mean(), factors.mean()))[:, None] - rf
    )
    p, q = np.vstack(rows), np.vstack(targets)
    c = sa / len(assets)
    noise = np.diag(np.diag(p @ c @ p.T))
    gain = np.linalg.solve(p @ c @ p.T + noise, p @ c).T
    posterior = prior + gain @ (q - p @ prior)
    posterior_cov = sa + c - gain @ p @ c
    expected_weights = np.linalg.solve(delta * posterior_cov, posterior)[:3]
    np.testing.assert_allclose(mu, (posterior[:3] + rf).T, rtol=1e-9, atol=1e-12)
    np.testing.assert_allclose(cov, posterior_cov[:3, :3], rtol=1e-9, atol=1e-12)
    np.testing.assert_allclose(weights, expected_weights, rtol=1e-9, atol=1e-12)
