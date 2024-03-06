""""""  #
"""
Copyright (c) 2020-2022, Dany Cajas
All rights reserved.
This work is licensed under BSD 3-Clause "New" or "Revised" License.
License available at https://github.com/dcajasn/Riskfolio-Lib/blob/master/LICENSE.txt
"""

import os
import numpy as np
import pandas as pd
import riskfolio as rp

assets = ["JCI", "TGT", "CMCSA", "CPB", "MO", "AMZN", "APA", "MMC", "JPM", "ZION"]
assets.sort()
benchmark = ["SPY"]


def resource(name):
    return os.path.join(os.path.abspath(os.path.dirname(__file__)), name)


def get_data(name):
    return pd.read_csv(resource(name), parse_dates=True, index_col=0)


def test_classic_minrisk_optimization():

    Y = get_data("stock_prices.csv")
    Y = Y[assets].pct_change().dropna().iloc[-200:]

    port = rp.Portfolio(returns=Y)

    method_mu = "hist"
    method_cov = "hist"

    port.assets_stats(method_mu=method_mu, method_cov=method_cov, d=0.94)
    port.alpha = 0.05
    port.solvers = ['CLARABEL', 'SCS', 'ECOS']

    model = "Classic"
    obj = "MinRisk"
    hist = True
    rf = 0
    l = 0

    rms = [
        "MV",
        "MAD",
        "GMD",
        "MSV",
        "FLPM",
        "SLPM",
        "CVaR",
        "TG",
        "EVaR",
        "WR",
        "RG",
        "CVRG",
        "TGRG",
        "MDD",
        "ADD",
        "CDaR",
        "EDaR",
        "UCI",
    ]

    w_1 = pd.DataFrame([])

    for i in rms:
        w = port.optimization(model=model, rm=i, obj=obj, rf=rf, l=l, hist=hist)
        w_1 = pd.concat([w_1, w], axis=1)

    w_1.columns = rms
    # w_1.to_csv("Classic_MinRisk.csv")

    w_2 = get_data("Classic_MinRisk.csv")

    a = np.testing.assert_array_almost_equal(w_1.to_numpy(), w_2.to_numpy(), decimal=6)
    if a is None:
        print("There are no errors in test_classic_minrisk_optimization")


def test_classic_sharpe_optimization():

    Y = get_data("stock_prices.csv")
    Y = Y[assets].pct_change().dropna().iloc[-200:]

    port = rp.Portfolio(returns=Y)

    method_mu = "hist"
    method_cov = "hist"

    port.assets_stats(method_mu=method_mu, method_cov=method_cov, d=0.94)
    port.alpha = 0.05
    port.solvers = ['CLARABEL', 'SCS', 'ECOS']

    model = "Classic"
    obj = "Sharpe"
    hist = True
    rf = 0
    l = 0

    rms = [
        "MV",
        "MAD",
        "GMD",
        "MSV",
        "FLPM",
        "SLPM",
        "CVaR",
        "TG",
        "EVaR",
        "WR",
        "RG",
        "CVRG",
        "TGRG",
        "MDD",
        "ADD",
        "CDaR",
        "EDaR",
        "UCI",
    ]

    w_1 = pd.DataFrame([])

    for i in rms:
        w = port.optimization(model=model, rm=i, obj=obj, rf=rf, l=l, hist=hist)
        w_1 = pd.concat([w_1, w], axis=1)

    w_1.columns = rms
    # w_1.to_csv("Classic_Sharpe.csv")

    w_2 = get_data("Classic_Sharpe.csv")

    a = np.testing.assert_array_almost_equal(w_1.to_numpy(), w_2.to_numpy(), decimal=6)
    if a is None:
        print("There are no errors in test_classic_sharpe_optimization")


def test_classic_riskparity_optimization():

    Y = get_data("stock_prices.csv")
    Y = Y[assets].pct_change().dropna().iloc[-200:]

    port = rp.Portfolio(returns=Y)

    method_mu = "hist"
    method_cov = "hist"

    port.assets_stats(method_mu=method_mu, method_cov=method_cov, d=0.94)
    port.alpha = 0.05
    port.solvers = ['CLARABEL', 'ECOS', 'SCS']

    model = "Classic"
    hist = True
    rf = 0
    b = None

    rms = [
        "MV",
        "MAD",
        "GMD",
        "MSV",
        "FLPM",
        "SLPM",
        "CVaR",
        "TG",
        "CVRG",
        "TGRG",
        "EVaR",
        "CDaR",
        "EDaR",
        "UCI",
    ]

    w_1 = pd.DataFrame([])

    for i in rms:
        w = port.rp_optimization(model=model, rm=i, rf=rf, b=b, hist=hist)
        w_1 = pd.concat([w_1, w], axis=1)

    w_1.columns = rms
    # w_1.to_csv("Classic_RP.csv")

    w_2 = get_data("Classic_RP.csv")

    a = np.testing.assert_array_almost_equal(w_1.to_numpy(), w_2.to_numpy(), decimal=6)
    if a is None:
        print("There are no errors in test_classic_riskparity_optimization")


def test_worst_case_optimization():

    Y = get_data("stock_prices.csv")
    Y = Y[assets].pct_change().dropna().iloc[-200:]

    port = rp.Portfolio(returns=Y)

    method_mu = "hist"
    method_cov = "hist"

    port.assets_stats(method_mu=method_mu, method_cov=method_cov, d=0.94)
    port.solvers = ['CLARABEL', 'ECOS', 'SCS']

    box = 's'
    ellip = 's'
    q = 0.05
    n_sim = 3000
    window = 3
    dmu = 0.1
    dcov = 0.1
    seed = 0

    port.wc_stats(box=box, ellip=ellip, q=q, n_sim=n_sim, window=window, dmu=dmu, dcov=dcov, seed=seed)

    Umus = ['box', 'ellip']
    Ucovs = ['box', 'ellip']
    objs = ['MinRisk', 'Sharpe']
    rf = 0
    l = 0

    w_1 = pd.DataFrame([])
    headers = []
    for obj in objs:
        for Umu in Umus:
            for Ucov in Ucovs:
                w = port.wc_optimization(obj=obj, rf=rf, l=l, Umu=Umu, Ucov=Ucov)
                w_1 = pd.concat([w_1, w], axis=1)
                headers += [obj + '-' + Umu + '-' + Ucov]

    w_1.columns = headers
    # w_1.to_csv("Classic_WC.csv")

    w_2 = get_data("Classic_WC.csv")

    a = np.testing.assert_array_almost_equal(w_1.to_numpy(), w_2.to_numpy(), decimal=6)
    if a is None:
        print("There are no errors in test_worst_case_optimization")


def test_hc_hrp_optimization():

    Y = get_data("stock_prices.csv")
    Y = Y[assets].pct_change().dropna().iloc[-200:]

    port = rp.HCPortfolio(returns=Y)

    model = "HRP"
    codependence = "pearson"
    rf = 0
    linkage = "single"
    max_k = 10
    leaf_order = True

    rms = [
        "vol",
        "MV",
        "MAD",
        "GMD",
        "MSV",
        "FLPM",
        "SLPM",
        "VaR",
        "CVaR",
        "TG",
        "EVaR",
        "WR",
        "RG",
        "CVRG",
        "TGRG",
        "MDD",
        "ADD",
        "DaR",
        "CDaR",
        "EDaR",
        "UCI",
        "MDD_Rel",
        "ADD_Rel",
        "DaR_Rel",
        "CDaR_Rel",
        "EDaR_Rel",
        "UCI_Rel",
    ]

    w_1 = pd.DataFrame([])

    for i in rms:
        w = port.optimization(
            model=model,
            codependence=codependence,
            rm=i,
            rf=rf,
            linkage=linkage,
            max_k=max_k,
            leaf_order=leaf_order,
        )

        w_1 = pd.concat([w_1, w], axis=1)

    w_1.columns = rms

    w_2 = get_data("HC_HRP.csv")

    a = np.testing.assert_array_almost_equal(w_1.to_numpy(), w_2.to_numpy(), decimal=6)
    if a is None:
        print("There are no errors in test_hc_hrp_optimization")


def test_hc_herc_optimization():

    Y = get_data("stock_prices.csv")
    Y = Y[assets].pct_change().dropna().iloc[-200:]

    port = rp.HCPortfolio(returns=Y)

    model = "HERC"
    codependence = "pearson"
    rf = 0
    linkage = "ward"
    max_k = 10
    leaf_order = True

    rms = [
        "vol",
        "MV",
        "MAD",
        "GMD",
        "MSV",
        "FLPM",
        "SLPM",
        "VaR",
        "CVaR",
        "TG",
        "EVaR",
        "WR",
        "RG",
        "CVRG",
        "TGRG",
        "MDD",
        "ADD",
        "DaR",
        "CDaR",
        "EDaR",
        "UCI",
        "MDD_Rel",
        "ADD_Rel",
        "DaR_Rel",
        "CDaR_Rel",
        "EDaR_Rel",
        "UCI_Rel",
    ]

    w_1 = pd.DataFrame([])

    for i in rms:
        w = port.optimization(
            model=model,
            codependence=codependence,
            rm=i,
            rf=rf,
            linkage=linkage,
            max_k=max_k,
            leaf_order=leaf_order,
        )

        w_1 = pd.concat([w_1, w], axis=1)

    w_1.columns = rms

    w_2 = get_data("HC_HERC.csv")

    a = np.testing.assert_array_almost_equal(w_1.to_numpy(), w_2.to_numpy(), decimal=6)
    if a is None:
        print("There are no errors in test_hc_herc_optimization")


def test_hc_nco_optimization():

    Y = get_data("stock_prices.csv")
    Y = Y[assets].pct_change().dropna().iloc[-200:]

    port = rp.HCPortfolio(returns=Y)

    model = "NCO"
    codependence = "pearson"
    covariance = "hist"
    obj = "MinRisk"
    rf = 0
    linkage = "ward"
    max_k = 10
    leaf_order = True

    rms = [
        "MV",
        "MAD",
        "MSV",
        "FLPM",
        "SLPM",
        "CVaR",
        "EVaR",
        "WR",
        "MDD",
        "ADD",
        "CDaR",
        "EDaR",
        "UCI",
    ]

    w_1 = pd.DataFrame([])

    for i in rms:
        w = port.optimization(
            model=model,
            codependence=codependence,
            covariance=covariance,
            obj=obj,
            rm=i,
            rf=rf,
            linkage=linkage,
            max_k=max_k,
            leaf_order=leaf_order,
        )

        w_1 = pd.concat([w_1, w], axis=1)

    w_1.columns = rms

    w_2 = get_data("HC_NCO.csv")

    a = np.testing.assert_array_almost_equal(w_1.to_numpy(), w_2.to_numpy(), decimal=6)
    if a is None:
        print("There are no errors in test_hc_nco_optimization")


if __name__ == '__main__':
    test_classic_minrisk_optimization()
    test_classic_sharpe_optimization()
    test_classic_riskparity_optimization()
    test_worst_case_optimization()
    test_hc_hrp_optimization()
    test_hc_herc_optimization()
    test_hc_nco_optimization()