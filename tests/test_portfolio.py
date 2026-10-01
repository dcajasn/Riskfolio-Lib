""""""  #

# Copyright (c) 2020-2026, Dany Cajas
# All rights reserved.
# This work is licensed under BSD 3-Clause "New" or "Revised" License.
# License available at https://github.com/dcajasn/Riskfolio-Lib/blob/master/LICENSE.txt

import os

import numpy as np
import pandas as pd
import pytest

import riskfolio as rp

assets = ["JCI", "TGT", "CMCSA", "CPB", "MO", "AMZN", "APA", "MMC", "JPM", "ZION"]
assets.sort()
benchmark = ["SPY"]

# The recorded weights come from convex solvers, so they are not bit stable
# across operating systems, BLAS builds and solver releases. The previous suite
# asserted decimal=6, which no solver based test can hold.
#
# 1e-3 was chosen by measurement rather than taste. Swapping CLARABEL for SCS,
# the largest move that is *not* a regression is:
#
#   Classic_*   5.4e-4   several of these problems have many optima, so the
#                        weights move while the objective is unchanged to 1e-7
#                        relative (MSV, CVaR, CVRG, ADD, MAD)
#   HC_HRP      2.9e-4   26 of the 27 columns are bit identical; only EVaR,
#   HC_HERC     6.7e-5   which calls a solver to evaluate the risk, moves
#
# and a real change is far larger: re-recording these files after four major
# releases moved them by 0.0052 (smallest) to 0.37 (largest). So 1e-3 sits
# about 2x above solver noise and 5x below a genuine regression.
#
# What this does not catch: a model error small enough to move the weights by
# less than 1e-3. Empirically, changing the inverse variance exponent in the
# HRP naive risk allocation from 2 to 1.9 is caught; 1.95 is not.
ATOL = 1e-3


def resource(name):
    return os.path.join(os.path.abspath(os.path.dirname(__file__)), name)


def get_data(name):
    return pd.read_csv(resource(name), parse_dates=True, index_col=0)


def returns():
    Y = get_data("stock_prices.csv")
    return Y[assets].pct_change().dropna().iloc[-200:]


def _portfolio(solvers):
    port = rp.Portfolio(returns=returns())
    port.assets_stats(
        method_mu="hist",
        method_cov="hist",
        dict_mu={"d": 0.94},
        dict_cov={"d": 0.94},
    )
    port.alpha = 0.05
    port.solvers = list(solvers)
    return port


def _concat(frames, columns):
    out = pd.concat(frames, axis=1)
    out.columns = columns
    return out


###############################################################################
# The computation behind each recorded file.
#
# The tests and tests/regenerate_goldens.py both call these, so a recorded file
# cannot drift away from the code that produced it.
###############################################################################


def classic_minrisk():
    port = _portfolio(["CLARABEL", "SCS", "ECOS"])
    rms = [
        "MV", "MAD", "GMD", "MSV", "FLPM", "SLPM", "CVaR", "TG", "EVaR",
        "WR", "RG", "CVRG", "TGRG", "MDD", "ADD", "CDaR", "EDaR", "UCI",
    ]
    ws = [
        port.optimization(model="Classic", rm=rm, obj="MinRisk", rf=0, l=0, hist=True)
        for rm in rms
    ]
    return _concat(ws, rms)


def classic_sharpe():
    port = _portfolio(["CLARABEL", "SCS", "ECOS"])
    rms = [
        "MV", "MAD", "GMD", "MSV", "FLPM", "SLPM", "CVaR", "TG", "EVaR",
        "WR", "RG", "CVRG", "TGRG", "MDD", "ADD", "CDaR", "EDaR", "UCI",
    ]
    ws = [
        port.optimization(model="Classic", rm=rm, obj="Sharpe", rf=0, l=0, hist=True)
        for rm in rms
    ]
    return _concat(ws, rms)


def classic_riskparity():
    port = _portfolio(["CLARABEL", "ECOS", "SCS"])
    rms = [
        "MV", "MAD", "GMD", "MSV", "FLPM", "SLPM", "CVaR", "TG", "CVRG",
        "TGRG", "EVaR", "CDaR", "EDaR", "UCI",
    ]
    ws = [
        port.rp_optimization(model="Classic", rm=rm, rf=0, b=None, hist=True)
        for rm in rms
    ]
    return _concat(ws, rms)


def classic_worst_case():
    port = _portfolio(["CLARABEL", "ECOS", "SCS"])
    port.wc_stats(
        box="s", ellip="s", q=0.05, n_sim=3000, window=3, dmu=0.1, dcov=0.1, seed=0
    )
    ws, headers = [], []
    for obj in ["MinRisk", "Sharpe"]:
        for Umu in ["box", "ellip"]:
            for Ucov in ["box", "ellip"]:
                ws.append(port.wc_optimization(obj=obj, rf=0, l=0, Umu=Umu, Ucov=Ucov))
                headers.append(f"{obj}-{Umu}-{Ucov}")
    return _concat(ws, headers)


_HC_RMS = [
    "vol", "MV", "MAD", "GMD", "MSV", "FLPM", "SLPM", "VaR", "CVaR", "TG",
    "EVaR", "WR", "RG", "CVRG", "TGRG", "MDD", "ADD", "DaR", "CDaR", "EDaR",
    "UCI", "MDD_Rel", "ADD_Rel", "DaR_Rel", "CDaR_Rel", "EDaR_Rel", "UCI_Rel",
]


def _hc(model, linkage, rms, **kwargs):
    port = rp.HCPortfolio(returns=returns())
    ws = [
        port.optimization(
            model=model,
            codependence="pearson",
            rm=rm,
            rf=0,
            linkage=linkage,
            max_k=10,
            leaf_order=True,
            **kwargs,
        )
        for rm in rms
    ]
    return _concat(ws, rms)


def hc_hrp():
    return _hc("HRP", "single", _HC_RMS)


def hc_herc():
    return _hc("HERC", "ward", _HC_RMS)


def hc_nco():
    rms = [
        "MV", "MAD", "MSV", "FLPM", "SLPM", "CVaR", "EVaR", "WR", "MDD",
        "ADD", "CDaR", "EDaR", "UCI",
    ]
    return _hc("NCO", "ward", rms, obj="MinRisk", method_cov="hist")


GOLDENS = {
    "Classic_MinRisk.csv": classic_minrisk,
    "Classic_Sharpe.csv": classic_sharpe,
    "Classic_RP.csv": classic_riskparity,
    "Classic_WC.csv": classic_worst_case,
    "HC_HRP.csv": hc_hrp,
    "HC_HERC.csv": hc_herc,
    "HC_NCO.csv": hc_nco,
}


def check_invariants(name, w):
    """Properties that must hold whatever the recorded numbers say."""
    assert list(w.index) == assets, f"{name}: unexpected asset order"
    assert np.isfinite(w.to_numpy()).all(), f"{name}: non finite weights"
    assert (w.to_numpy() >= -1e-6).all(), (
        f"{name}: negative weights in a long only model"
    )
    budget = w.sum(axis=0).to_numpy()
    np.testing.assert_allclose(
        budget, 1.0, atol=1e-4, err_msg=f"{name}: weights do not sum to one"
    )


@pytest.mark.parametrize("name", sorted(GOLDENS))
def test_optimization_matches_recorded_weights(name):
    w = GOLDENS[name]()
    check_invariants(name, w)

    expected = get_data(name)
    assert list(w.columns) == list(expected.columns), f"{name}: column set changed"
    assert list(w.index) == list(expected.index), f"{name}: asset set changed"
    np.testing.assert_allclose(
        w.to_numpy(),
        expected.to_numpy(),
        atol=ATOL,
        rtol=0,
        err_msg=(
            f"{name} differs from the recorded weights. If this is a deliberate "
            f"change, rerun tests/regenerate_goldens.py and review the diff."
        ),
    )
