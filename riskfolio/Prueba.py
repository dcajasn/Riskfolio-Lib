#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 00:13:42 2020

@author: danycajas
"""


import numpy as np
import pandas as pd
import yfinance as yf

# # Libraries final
# import riskfolio.ParamsEstimation as pe
# import riskfolio.Portfolio as pf
# import riskfolio.PlotFunctions as plf

# Libraries test
import ParamsEstimation as pe
import Portfolio as pf
import PlotFunctions as plf

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

loadings = pe.loadings_matrix(factors, assets)

#%%

X = factors
Y = assets

# Creando el Portafolio
port = pf.Portfolio(returns=Y)
pd.options.display.float_format = "{:.4%}".format

rm = "MAD"
method_mu = "hist"
method_cov = "hist"

# Prueba modelos clasicos (solo data de retornos de activos)
port.assets_stats(method_mu=method_mu, method_cov=method_cov, d=0.94)

w1 = port.optimization(model="Classic", rm=rm, obj="Sharpe", rf=0.0, l=0, hist=True)

# Prueba modelos de factores (data de retornos de activos y factores)
port.factors = X
port.factors_stats()

w2 = port.optimization(model="FM", rm="MV", obj="Sharpe", rf=0.0, l=3, hist=False)
print(w2)

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
    w=w1,
    label="Portfolio",
    marker="*",
    s=16,
    c="r",
    height=6,
    width=10,
    ax=None,
)

# ax = plf.plot_pie(w=w1, title='Portafolio', height=6, width=10, cmap="tab20",
#                   ax=None)

# # # w0 = pd.DataFrame(np.ones((100,1))/100)

# # # ax = plf.plot_pie(w=w0, title='Portafolio', height=6, width=10, cmap = "tab20", ax=None)

# # ax = plf.plot_frontier_area(w_frontier=w_5, cmap="tab20", height=6, width=10, ax=None)

# ax = plf.plot_hist(data=Y, w=w1, alpha=0.01, bins=50, height=6, width=10,
#                     ax=None)

# ax = plf.plot_drawdown(data=Y, w=w1, alpha=0.01, height=8, width=10, ax=None)

# # ax = plf.plot_series(data=Y, w=w_5, cmap='tab20', height=6, width=10, ax=None)
