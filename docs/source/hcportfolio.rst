##############################################
Hierarchical Clustering Portfolio Optimization
##############################################

.. raw:: html

    <a href='https://ko-fi.com/B0B833SXD' target='_blank'><img height='36'style='border:0px;height:36px;' src='https://cdn.ko-fi.com/cdn/kofi1.png?v=2' border='0' alt='Buy Me a Coffee at ko-fi.com' /></a>

Some Theory
===========

Hierarchical Clustering Portfolio Optimization
----------------------------------------------

Riskfolio-Lib allows to calculate the new machine learning asset allocation models. The available models are:

- Hierarchical Risk Parity (HRP) :cite:`c-Prado1`, :cite:`c-Sjostrand`.
- Hierarchical Equal Risk Contribution (HERC) :cite:`c-Raffinot2`, :cite:`c-Sjostrand`.

In both cases we have the option to use the following 22 risk measures to calculate HRP and HERC portfolios using naive risk parity:

    - Standard Deviation.
    - Variance.
    - Semi Standard Deviation.
    - Mean Absolute Deviation (MAD).
    - First Lower Partial Moment (Omega Ratio).
    - Second Lower Partial Moment (Sortino Ratio).
    - Value at Risk (VaR).
    - Conditional Value at Risk (CVaR).
    - Entropic Value at Risk (EVaR).
    - Worst Case Realization (Minimax Model).
    - Maximum Drawdown (Calmar Ratio) for compounded and uncompounded cumulative returns.
    - Average Drawdown for compounded and uncompounded cumulative returns.
    - Drawdown at Risk (DaR) for compounded and uncompounded cumulative returns.
    - Conditional Drawdown at Risk (CDaR) for compounded and uncompounded cumulative returns.
    - Entropic Drawdown at Risk (EDaR) for compounded and uncompounded cumulative returns.
    - Ulcer Index for compounded and uncompounded cumulative returns.


Module Methods
==============

.. automodule:: HCPortfolio
   :members:
   :private-members:


Bibliography
============

.. bibliography:: biblio.bib
   :style: unsrt
   :labelprefix: C
   :keyprefix: c-