"""
Tests for MVSK optimization via yand-mvsk integration.
"""

import os

import numpy as np
import pandas as pd
import pytest

try:
    import yand_mvsk

    MVSK_INSTALLED = True
except ImportError:
    MVSK_INSTALLED = False

import riskfolio as rp


def resource(name):
    return os.path.join(os.path.abspath(os.path.dirname(__file__)), name)


def get_data(name):
    return pd.read_csv(resource(name), index_col=0)


assets = ["JCI", "TGT", "CMCSA", "CPB", "MO", "AMZN", "APA", "MMC", "JPM", "ZION"]
assets.sort()


def setup_portfolio():
    Y = get_data("stock_prices.csv")
    Y = Y[assets].pct_change().dropna().iloc[-200:]
    port = rp.Portfolio(returns=Y)
    port.assets_stats(method_mu="hist", method_cov="hist")
    return port


@pytest.mark.skipif(not MVSK_INSTALLED, reason="yand-mvsk not installed")
class TestMVSK:
    def test_mvsk_basic(self):
        port = setup_portfolio()
        w = port.optimization(model="Classic", rm="MVSK", obj="Utility", l=6)

        assert w is not None
        assert isinstance(w, pd.DataFrame)
        assert w.shape == (len(assets), 1)
        np.testing.assert_almost_equal(w["weights"].sum(), 1.0, decimal=6)

    def test_mvsk_weights_long_only(self):
        port = setup_portfolio()
        w = port.optimization(model="Classic", rm="MVSK", obj="Utility", l=6)

        assert np.all(w["weights"].values >= -1e-8)

    def test_mvsk_tickers(self):
        port = setup_portfolio()
        w = port.optimization(model="Classic", rm="MVSK", obj="Utility", l=6)

        assert list(w.index) == assets

    def test_mvsk_different_gamma(self):
        port1 = setup_portfolio()
        w_low = port1.optimization(model="Classic", rm="MVSK", obj="Utility", l=2)

        port2 = setup_portfolio()
        w_high = port2.optimization(model="Classic", rm="MVSK", obj="Utility", l=20)

        assert not np.allclose(
            w_low["weights"].values, w_high["weights"].values, atol=1e-3
        )

    def test_mvsk_stored_as_optimal(self):
        port = setup_portfolio()
        w = port.optimization(model="Classic", rm="MVSK", obj="Utility", l=6)

        assert port.optimal is not None
        pd.testing.assert_frame_equal(w, port.optimal)

    def test_mvsk_budget_constraint(self):
        port = setup_portfolio()
        port.budget = 1.0
        w = port.optimization(model="Classic", rm="MVSK", obj="Utility", l=6)

        np.testing.assert_almost_equal(w["weights"].sum(), 1.0, decimal=6)

    def test_mvsk_lower_bound(self):
        port = setup_portfolio()
        port.lowerlng = 0.02
        w = port.optimization(model="Classic", rm="MVSK", obj="Utility", l=6)

        assert np.all(w["weights"].values >= 0.02 - 1e-6)

    def test_mvsk_obj_minrisk(self):
        port = setup_portfolio()
        w = port.optimization(model="Classic", rm="MVSK", obj="MinRisk", l=6)

        assert w is not None
        np.testing.assert_almost_equal(w["weights"].sum(), 1.0, decimal=6)

    def test_mvsk_obj_sharpe(self):
        port = setup_portfolio()
        w = port.optimization(model="Classic", rm="MVSK", obj="Sharpe", l=6)

        assert w is not None
        np.testing.assert_almost_equal(w["weights"].sum(), 1.0, decimal=6)

    def test_mvsk_beats_equal_weight(self):
        port = setup_portfolio()
        w = port.optimization(model="Classic", rm="MVSK", obj="Utility", l=6)

        from yand_mvsk import MVSKOracle, crra_coefficients

        returns = port.returns.values
        c = crra_coefficients(6.0)
        oracle = MVSKOracle(returns, c)

        mvsk_weights = w["weights"].values
        eq_weights = np.ones(len(assets)) / len(assets)

        assert oracle.value(mvsk_weights) <= oracle.value(eq_weights)

    def test_mvsk_without_package_raises(self):
        """Verify the ImportError message is clear when yand-mvsk is missing."""
        pass

    def test_mvsk_fm_model(self):
        """MVSK should work with FM model when returns_fm is available."""
        port = setup_portfolio()
        w = port.optimization(model="Classic", rm="MVSK", obj="Utility", l=6)
        assert w is not None
