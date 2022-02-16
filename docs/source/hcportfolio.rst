##############################################
Hierarchical Clustering Portfolio Optimization
##############################################

.. raw:: html

    <a href='https://ko-fi.com/B0B833SXD' target='_blank'><img height='36'style='border:0px;height:36px;' src='https://cdn.ko-fi.com/cdn/kofi1.png?v=2' border='0' alt='Buy Me a Coffee at ko-fi.com' /></a><br>

Some Theory
===========

Hierarchical Clustering Portfolio Optimization
----------------------------------------------

Riskfolio-Lib allows to calculate the new machine learning asset allocation models. The available models are:

- Hierarchical Risk Parity (HRP) :cite:`c-Prado1`, :cite:`c-Sjostrand`, :cite:`c-Pfitzinger`.
- Hierarchical Equal Risk Contribution (HERC) :cite:`c-Raffinot2`, :cite:`c-Sjostrand`.
- Nested Clustered Optimization (NCO) :cite:`c-Prado2`, :cite:`c-Sjostrand`.

In the first two cases we have the option to use the following 27 risk measures to calculate HRP and HERC portfolios using naive risk parity:

**Dispersion Risk Measures:**

- Standard Deviation.
- Variance.
- Mean Absolute Deviation (MAD).
- Gini Mean Difference (GMD).
- Range.
- Conditional Value at Risk Range.
- Tail Gini Range.

**Downside Risk Measures:**

- Semi Standard Deviation.
- First Lower Partial Moment (Omega Ratio).
- Second Lower Partial Moment (Sortino Ratio).
- Value at Risk (VaR).
- Conditional Value at Risk (CVaR).
- Entropic Value at Risk (EVaR).
- Tail Gini.
- Worst Case Realization (Minimax).

**Drawdown Risk Measures:**

- Maximum Drawdown (Calmar Ratio) for compounded and uncompounded cumulative returns.
- Average Drawdown for compounded and uncompounded cumulative returns.
- Drawdown at Risk (DaR) for compounded and uncompounded cumulative returns.
- Conditional Drawdown at Risk (CDaR) for compounded and uncompounded cumulative returns.
- Entropic Drawdown at Risk (EDaR) for compounded and uncompounded cumulative returns.
- Ulcer Index for compounded and uncompounded cumulative returns.

For the NCO model we have the option to use four objective functions with the available risk measures to each objective:

- Minimize the selected risk measure.
- Maximize the Utility function :math:`\mu w - l \phi_{i}(w)`.
- Maximize the risk adjusted return ratio based on the selected risk measure.
- Equally risk contribution portfolio of the selected risk measure.


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