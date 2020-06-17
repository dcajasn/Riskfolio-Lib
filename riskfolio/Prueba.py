#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 00:13:42 2020

@author: danycajas
"""


import numpy as np
import pandas as pd
import yfinance as yf

# Libraries final
import riskfolio.ParamsEstimation as pe
import riskfolio.Portfolio as pf
import riskfolio.PlotFunctions as plf

# # Libraries test
# import ParamsEstimation as pe
# import Portfolio as pf
# import PlotFunctions as plf

#%%

yf.pdr_override()

# Fecha de inicio y final
start = "2017-01-01"
end = "2019-12-30"

# Data para la prueba
# assets = ['FB', 'GOOGL', 'NFLX', 'BAC', 'WFC', 'TLT', 'SHV']
assets = [
    "JCI",
    "TGT",
    "CMCSA",
    "CPB",
    "MO",
    "NBL",
    "APA",
    "MMC",
    "JPM",
    "ZION",
    "PSA",
    "AGN",
    "BAX",
    "BMY",
    "LUV",
    "PCAR",
    "TXT",
    "DHR",
    "DE",
    "MSFT",
    "HPQ",
    "SEE",
    "VZ",
    "CNP",
    "NI",
]
assets.sort()

factors = ["MTUM", "QUAL", "VLUE", "SIZE", "USMV"]
factors.sort()

tickers = assets + factors
tickers.sort()

# Descargando la data de precios de activos y factores
data = yf.download(tickers + factors, start=start, end=end)
data = data.loc[:, ("Adj Close", slice(None))]
data.columns = tickers

assets = data[assets].pct_change().dropna()
factors = data[factors].pct_change().dropna()

loadings = pe.loadings_matrix(
    factors, assets, feature_selection="PCR", n_components=0.90
)

mu, cov, returns, nav = pe.risk_factors(
    factors,
    assets,
    B=None,
    method_mu="hist",
    method_cov="hist",
    threshold=0.05,
    feature_selection="PCR",
    stepwise="Forward",
    n_components=0.95,
    error=True,
)

#%%

X = factors
Y = assets

# Creando el Portafolio
port = pf.Portfolio(returns=Y)

# Pesos en Corto
port.sht = True
port.uppersht = 0.2

pd.options.display.float_format = "{:.4%}".format

rm = "MAD"
method_mu = "hist"
method_cov = "hist"

# Prueba modelos clasicos (solo data de retornos de activos)
port.assets_stats(method_mu=method_mu, method_cov=method_cov, d=0.94)

# w1 = port.optimization(model="Classic", rm=rm, obj="Sharpe", rf=0.0, l=0, hist=True)

# # Prueba modelos de factores (data de retornos de activos y factores)
# port.factors = X
# port.factors_stats()

# w2 = port.optimization(model="FM", rm="MV", obj="Sharpe", rf=0.0, l=3, hist=False)
# print(w2)

# Prueba modelos de factores (data de retornos de activos y factores)
port.factors = X
port.factors_stats(feature_selection="stepwise", n_components=0.90)

w2_1 = port.optimization(model="FM", rm="MV", obj="Sharpe", rf=0.0, l=3, hist=False)
print(w2_1)


# # # Prueba modelo de Black Litterman
# # n_views = 3
# # n_assets = w1.shape[0]

# # delta = 2
# # P = np.zeros((n_views, n_assets))
# # Q = np.zeros((n_views, 1))

# # P[0, 1] = 1
# # P[0, 5] = -1
# # P[1, 3] = 1
# # P[1, 9] = -1
# # P[2, 0] = 1

# # Q[0, 0] = 0.05/252
# # Q[1, 0] = 0.09/252
# # Q[2, 0] = 0.07/252

# # port.blacklitterman_stats(P, Q, delta, rf=0, w=w1, equilibrium=False)

# # w3 = port.optimization(model='BL', rm='MV', obj='Sharpe', rf=0.0, l=3, hist=False)
# # print(w3)
# # ws = pd.concat([w1, w2, w3], axis=1)

# # # Prueba modelos de Black Litterman con Factores


# # Prueba Limites Frontera Eficiente
# w_4 = port.frontier_limits(model='Classic', rm=rm, rf=0,
#                            hist=True)


# Prueba Frontera Efficiente

w_5 = port.efficient_frontier(model="Classic", rm=rm, points=20, rf=0, hist=True)


#%%

label = "Max Risk Adjusted Return Portfolio"
mu = port.mu
cov = port.cov
returns = port.returns

ax = plf.plot_frontier(
    w_frontier=w_5,
    mu=mu,
    cov=cov,
    returns=returns,
    rm=rm,
    rf=0,
    alpha=0.01,
    cmap="viridis",
    w=w2_1,
    label="Portfolio",
    marker="*",
    s=16,
    c="r",
    height=6,
    width=10,
    ax=None,
)

ax = plf.plot_pie(w=w2_1, title="Portafolio", height=6, width=10, cmap="tab20", ax=None)

# # # w0 = pd.DataFrame(np.ones((100,1))/100)

# # # ax = plf.plot_pie(w=w0, title='Portafolio', height=6, width=10, cmap = "tab20", ax=None)

# # ax = plf.plot_frontier_area(w_frontier=w_5, cmap="tab20", height=6, width=10, ax=None)

# ax = plf.plot_hist(data=Y, w=w1, alpha=0.01, bins=50, height=6, width=10,
#                     ax=None)

# ax = plf.plot_drawdown(data=Y, w=w1, alpha=0.01, height=8, width=10, ax=None)

# # ax = plf.plot_series(data=Y, w=w_5, cmap='tab20', height=6, width=10, ax=None)


#%%
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import statsmodels.api as sm

y = Y["BAX"]

criterion = "R2_A"
threshold = 0.05
verbose = True
included = []
value = 0

# while value <= threshold:
#     excluded = list(set(X.columns) - set(included))
#     best_pvalue = 999999
#     new_feature = None
#     for i in excluded:
#         factors = included + [i]
#         X1 = X[factors]
#         X1 = sm.add_constant(X1)
#         results = sm.OLS(y, X1).fit()
#         new_pvalues = results.pvalues
#         cond_1 = new_pvalues[new_pvalues.index != "const"].max()
#         if best_pvalue > new_pvalues[i] and cond_1 <= threshold:
#             best_pvalue = results.pvalues[i]
#             new_feature = i
#             pvalues = new_pvalues.copy()

#     value = pvalues[pvalues.index != "const"].max()

#     if new_feature is None:
#         break

#     included.append(new_feature)

#     if verbose:
#         print("Add {} with p-value {:.6}".format(new_feature, best_pvalue))


aic = 1e10
sic = 1e10
r2 = -1e10
r2_a = -1e10

excluded = X.columns.tolist()

for i in range(X.shape[1]):
    j = 0
    value = None
    for i in excluded:
        factors = included.copy()
        factors.append(i)
        X1 = X[factors]
        X1 = sm.add_constant(X1)
        results = sm.OLS(y, X1).fit()

        if criterion == "AIC":
            if results.aic < aic:
                value = i
                aic = results.aic
        if criterion == "SIC":
            if results.bic < sic:
                value = i
                sic = results.bic
        if criterion == "R2":
            if results.rsquared > r2:
                value = i
                r2 = results.rsquared
        if criterion == "R2_A":
            if results.rsquared_adj > r2_a:
                value = i
                r2_a = results.rsquared_adj

        j += 1
        if j == len(excluded):
            if value is None:
                break
            else:
                excluded.remove(value)
                included.append(value)
                if verbose:
                    if criterion == "AIC":
                        print("Add {} with AIC {:.6}".format(value, results.aic))
                    elif criterion == "SIC":
                        print("Add {} with SIC {:.6}".format(value, results.bic))
                    elif criterion == "R2":
                        print("Add {} with R2 {:.6}".format(value, results.rsquared))
                    elif criterion == "R2_A":
                        print(
                            "Add {} with Adjusted R2 {:.6}".format(
                                value, results.rsquared_adj
                            )
                        )
