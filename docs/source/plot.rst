.. Riskfolio documentation master file, created by
   sphinx-quickstart on Fri Feb  7 02:11:42 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

##############
Plot Functions
##############

.. raw:: html

    <a href='https://ko-fi.com/B0B833SXD' target='_blank'><img height='36'style='border:0px;height:36px;' src='https://cdn.ko-fi.com/cdn/kofi1.png?v=2' border='0' alt='Buy Me a Coffee at ko-fi.com' /></a>

This module has functions that allows us to create charts that 
helps us to analyze quickly the properties of our optimal portfolios.

The following example construct the portfolios and the efficient frontier that
will be plot using the functions of this module.

Example
=======

::
    
    import numpy as np
    import pandas as pd
    import yfinance as yf
    import riskfolio as rp
        
    # Date range
    start = '2016-01-01'
    end = '2019-12-30'

    # Tickers of assets
    tickers = ['JCI', 'TGT', 'CMCSA', 'CPB', 'MO', 'APA', 'MMC', 'JPM',
               'ZION', 'PSA', 'BAX', 'BMY', 'LUV', 'PCAR', 'TXT', 'TMO',
               'DE', 'MSFT', 'HPQ', 'SEE', 'VZ', 'CNP', 'NI', 'T', 'BA']
    tickers.sort()
    
    # Downloading the data
    data = yf.download(tickers, start = start, end = end)
    data = data.loc[:,('Adj Close', slice(None))]
    data.columns = tickers
    assets = data.pct_change().dropna()

    Y = assets
    
    # Creating the Portfolio Object
    port = rp.Portfolio(returns=Y)
    
    # To display dataframes values in percentage format
    pd.options.display.float_format = '{:.4%}'.format
    
    # Choose the risk measure
    rm = 'MSV'  # Semi Standard Deviation
    
    # Estimate inputs of the model (historical estimates)
    method_mu='hist' # Method to estimate expected returns based on historical data.
    method_cov='hist' # Method to estimate covariance matrix based on historical data.

    port.assets_stats(method_mu=method_mu, method_cov=method_cov, d=0.94)
    
    # Estimate the portfolio that maximizes the risk adjusted return ratio
    w1 = port.optimization(model='Classic', rm=rm, obj='Sharpe', rf=0.0, l=0, hist=True)

    # Estimate points in the efficient frontier mean - semi standard deviation
    ws = port.efficient_frontier(model='Classic', rm=rm, points=20, rf=0, hist=True)

    # Estimate the risk parity portfolio for semi standard deviation
    w2 = port.rp_optimization(model='Classic', rm=rm, rf=0, b=None, hist=True)



Module Functions
================

.. automodule:: PlotFunctions
   :members:
   :private-members:
