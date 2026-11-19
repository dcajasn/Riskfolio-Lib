import numpy as np
import pandas as pd

import riskfolio as rp
import riskfolio.src.AuxFunctions as af
import riskfolio.src.RiskFunctions as rk


def test_is_pos_def_returns_python_bool():
    assert af.is_pos_def(np.eye(3)) is True
    assert af.is_pos_def(np.array([[1.0, 2.0], [2.0, 1.0]])) is False


def test_hrp_constraints_applies_enabled_constraints():
    asset_classes = pd.DataFrame({"Assets": ["A", "B"]})
    constraints = pd.DataFrame(
        {
            "Disabled": [False],
            "Type": ["Assets"],
            "Set": [""],
            "Position": ["A"],
            "Sign": ["<="],
            "Weight": [0.3],
        }
    )
    w_max, _ = rp.hrp_constraints(constraints, asset_classes)
    assert w_max["A"] == 0.3


def _evar_hist_with_success(monkeypatch, success):
    minimize = rk.minimize
    found = {}

    def fake_minimize(*args, **kwargs):
        result = minimize(*args, **kwargs)
        found["z"] = result.x.item()
        result.success = success
        if not success:
            result.x = np.array([np.nan])
        return result

    monkeypatch.setattr(rk, "minimize", fake_minimize)
    X = np.random.default_rng(1).normal(0, 0.02, 500)
    value, z = rk.EVaR_Hist(X)
    return X, value, z, found["z"]


def test_evar_hist_uses_scipy_result_on_numpy_bool_success(monkeypatch):
    _, value, z, z_scipy = _evar_hist_with_success(monkeypatch, np.bool_(True))
    assert z == z_scipy
    assert np.isfinite(value)


def test_evar_hist_falls_back_to_cvxpy_on_numpy_bool_failure(monkeypatch):
    X, value, z, z_scipy = _evar_hist_with_success(monkeypatch, np.bool_(False))
    assert np.isfinite(z)
    assert np.isfinite(value)
    assert np.isclose(value, rk.Entropic_RM(X, z=z_scipy), rtol=1e-3)
