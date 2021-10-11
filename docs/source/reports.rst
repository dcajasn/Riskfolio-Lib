.. Riskfolio documentation master file, created by
   sphinx-quickstart on Fri Feb  7 02:11:42 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

#################
Reports Functions
#################

.. raw:: html

    <a href='https://ko-fi.com/B0B833SXD' target='_blank'><img height='36'style='border:0px;height:36px;' src='https://cdn.ko-fi.com/cdn/kofi1.png?v=2' border='0' alt='Buy Me a Coffee at ko-fi.com' /></a>

This section explains some functions that allows us to create Jupyter Notebook and Excel reports that helps us to analyze quickly the properties of our portfolios.

The following example build an optimum portfolio and create a Jupyter Notebook and Excel report using the functions of this module.

Example
=======

::
    
    import numpy as np
    import pandas as pd
    import yfinance as yf
    import riskfolio as rp
    
    yf.pdr_override()
    
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
    rm = 'MV'  # Standard Deviation
    
    # Estimate inputs of the model (historical estimates)
    method_mu='hist' # Method to estimate expected returns based on historical data.
    method_cov='hist' # Method to estimate covariance matrix based on historical data.

    port.assets_stats(method_mu=method_mu, method_cov=method_cov, d=0.94)
    
    # Estimate the portfolio that maximizes the risk adjusted return ratio
    w = port.optimization(model='Classic', rm=rm, obj='Sharpe', rf=0.0, l=0, hist=True)


Module Functions
================

.. automodule:: Reports
   :members:
   :private-members:
