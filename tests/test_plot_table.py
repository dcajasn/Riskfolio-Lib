import numpy as np
import pandas as pd

import riskfolio as rp
import riskfolio.src.RiskFunctions as rk


def test_plot_table_rows():
    index = pd.bdate_range("2024-01-01", periods=250)
    returns = pd.DataFrame(np.random.default_rng(0).normal(0, 0.01, (250, 2)), index)
    w = pd.DataFrame([0.6, 0.4])
    X = returns.to_numpy() @ w.to_numpy()

    cells = rp.plot_table(returns, w, alpha=0.1).tables[0].get_celld()
    rows = {
        cells[i, 0].get_text().get_text(): cells[i, 1].get_text().get_text()
        for i in range(1, 33)
    }

    tg = rk.TG(X, alpha=0.1) * 252**0.5
    assert rows["Tail Gini of Losses (TG) (2)"] == "{:.4%}".format(tg)
    assert rows["Drawdown at Risk (DaR)"] == "{:.4%}".format(rk.DaR_Abs(X, alpha=0.1))
