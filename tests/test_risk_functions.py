""""""  #

"""
Copyright (c) 2020-2026, Dany Cajas
All rights reserved.
This work is licensed under BSD 3-Clause "New" or "Revised" License.
License available at https://github.com/dcajasn/Riskfolio-Lib/blob/master/LICENSE.txt
"""

import numpy as np
import pytest
from scipy.stats import norm, skew, kurtosis

import riskfolio.src.RiskFunctions as rk


def _reference_mvar(x, alpha=0.05):
    """Independent Cornish-Fisher modified VaR, used to check MVaR_Hist."""
    x = np.asarray(x, dtype=float).flatten()
    mu = np.mean(x)
    sigma = np.std(x, ddof=1)
    S = skew(x, bias=False)
    K = kurtosis(x, fisher=True, bias=False)
    z = norm.ppf(alpha)
    z_cf = (
        z
        + (z**2 - 1) / 6 * S
        + (z**3 - 3 * z) / 24 * K
        - (2 * z**3 - 5 * z) / 36 * S**2
    )
    return -(mu + z_cf * sigma)


def test_mvar_hist_is_exported():
    assert "MVaR_Hist" in rk.__all__


def test_mvar_hist_matches_reference():
    rng = np.random.default_rng(123)
    for data in [rng.normal(0, 0.02, 500), rng.standard_t(5, 500) * 0.02 - 1e-3]:
        assert np.isclose(
            rk.MVaR_Hist(data, 0.05), _reference_mvar(data, 0.05), atol=1e-12
        )


def test_mvar_hist_reduces_to_gaussian_for_normal_data():
    rng = np.random.default_rng(0)
    x = rng.normal(0.001, 0.02, 50000)
    gaussian_var = -(np.mean(x) + norm.ppf(0.05) * np.std(x, ddof=1))
    assert np.isclose(rk.MVaR_Hist(x, 0.05), gaussian_var, atol=1e-3)


def test_mvar_hist_exceeds_gaussian_for_heavy_tails():
    # Symmetric heavy-tailed data (skew ~ 0, large positive excess kurtosis):
    # at the 1% level the Cornish-Fisher kurtosis term inflates the modified
    # VaR above the Gaussian VaR.
    rng = np.random.default_rng(7)
    x = rng.standard_t(5, 50000) * 0.01
    gaussian_var = -(np.mean(x) + norm.ppf(0.01) * np.std(x, ddof=1))
    assert rk.MVaR_Hist(x, 0.01) > gaussian_var


def test_mvar_hist_shape_handling():
    rng = np.random.default_rng(1)
    x = rng.normal(0, 0.02, 300)
    assert np.isclose(rk.MVaR_Hist(x.reshape(1, -1), 0.05), rk.MVaR_Hist(x, 0.05))
    with pytest.raises(ValueError):
        rk.MVaR_Hist(np.ones((3, 3)))


def test_mvar_hist_degenerate_inputs_are_finite():
    # Cornish-Fisher skew/kurtosis are undefined at zero variance and for a
    # single observation. MVaR_Hist should short-circuit to -mean (consistent
    # with VaR_Hist) and return a finite value rather than NaN.
    constant = np.full(100, 0.01)
    out = rk.MVaR_Hist(constant, 0.05)
    assert np.isfinite(out)
    assert np.isclose(out, -0.01)

    single = np.array([0.02])
    out_single = rk.MVaR_Hist(single, 0.05)
    assert np.isfinite(out_single)
    assert np.isclose(out_single, -0.02)
