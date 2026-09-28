import numpy as np
import pytest

import riskfolio.src.RiskFunctions as rk

RNG = np.random.default_rng(0)
CLEAN = RNG.normal(0.0005, 0.01, (250, 1))

# (function, extra args) for every single-series risk measure
FUNCS = [
    (rk.MAD, ()), (rk.SemiDeviation, ()), (rk.Kurtosis, ()), (rk.SemiKurtosis, ()),
    (rk.EvenMoment, ()), (rk.EvenSemiMoment, ()), (rk.VaR_Hist, ()), (rk.CVaR_Hist, ()),
    (rk.WR, ()), (rk.LPM, ()), (rk.Entropic_RM, ()), (rk.EVaR_Hist, ()),
    (rk.RLVaR_Hist, ()), (rk.MDD_Abs, ()), (rk.ADD_Abs, ()), (rk.DaR_Abs, ()),
    (rk.CDaR_Abs, ()), (rk.EDaR_Abs, ()), (rk.RLDaR_Abs, ()), (rk.UCI_Abs, ()),
    (rk.MDD_Rel, ()), (rk.ADD_Rel, ()), (rk.DaR_Rel, ()), (rk.CDaR_Rel, ()),
    (rk.EDaR_Rel, ()), (rk.RLDaR_Rel, ()), (rk.UCI_Rel, ()), (rk.GMD, ()),
    (rk.TG, ()), (rk.RG, ()), (rk.VRG, ()), (rk.CVRG, ()), (rk.TGRG, ()),
    (rk.EVRG, ()), (rk.RVRG, ()), (rk.L_Moment, ()), (rk.L_Moment_CRM, ()),
]


@pytest.mark.parametrize("func,args", FUNCS, ids=lambda f: getattr(f, "__name__", ""))
@pytest.mark.parametrize("bad", [np.nan, np.inf])
def test_non_finite_returns_are_rejected(func, args, bad):
    X = CLEAN.copy()
    X[100, 0] = bad
    with pytest.raises(ValueError, match="NaN or infinite"):
        func(X, *args)


def test_nan_no_longer_truncates_max_drawdown():
    # Before the guard, NAV turned NaN at the first missing value and every
    # later comparison was False, so the drawdown after it was ignored.
    X = np.array([[0.01], [-0.02], [np.nan], [0.03], [-0.05], [0.01]])
    with pytest.raises(ValueError):
        rk.MDD_Rel(X)
    assert rk.MDD_Rel(X[~np.isnan(X)].reshape(-1, 1)) == pytest.approx(0.05)


@pytest.mark.parametrize("func,args", FUNCS[:10], ids=lambda f: getattr(f, "__name__", ""))
def test_clean_returns_unaffected(func, args):
    assert np.isfinite(func(CLEAN, *args))
