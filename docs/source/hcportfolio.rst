##############################################
Hierarchical Clustering Portfolio Optimization
##############################################

.. raw:: html

    <a href="https://www.kqzyfj.com/click-101359873-15150084?url=https%3A%2F%2Flink.springer.com%2Fbook%2F9783031843037" target="_blank">
        <button style="padding:10px 20px; font-size:16px; background-color: #FFA500; color:white; border:none; border-radius:5px; cursor:pointer;">
            Buy Advanced Portfolio Optimization Book on Springer
        </button>
    </a>
    <br>
    <br>

.. image:: https://img.shields.io/static/v1?label=Sponsor&message=%E2%9D%A4&logo=GitHub&color=%23fe8e86
 :target: https://github.com/sponsors/dcajasn

.. raw:: html
   
    <br>
   
.. raw:: html

    <a href='https://ko-fi.com/B0B833SXD' target='_blank'><img height='36'style='border:0px;height:36px;' src='https://cdn.ko-fi.com/cdn/kofi1.png?v=2' border='0' alt='Buy Me a Coffee at ko-fi.com' /></a>

Some Theory
===========

Hierarchical Clustering Portfolio Optimization
----------------------------------------------

Riskfolio-Lib allows to calculate the new machine learning asset allocation models. The available models are:

- Hierarchical Risk Parity (HRP) :cite:`c-Prado1`, :cite:`c-Sjostrand`, :cite:`c-Pfitzinger`.
- Hierarchical Equal Risk Contribution (HERC) :cite:`c-Raffinot1`, :cite:`c-Raffinot2`, :cite:`c-Sjostrand`.
- Nested Clustered Optimization (NCO) :cite:`c-Prado2`, :cite:`c-Sjostrand`.

In the first two cases we have the option to use the following 32 risk measures to calculate HRP and HERC portfolios using naive risk parity:

**Dispersion Risk Measures:**

- Standard Deviation.
- Variance.
- Square Root Kurtosis.
- Mean Absolute Deviation (MAD).
- Gini Mean Difference (GMD).
- Conditional Value at Risk Range.
- Tail Gini Range.
- Range.

**Downside Risk Measures:**

- Semi Standard Deviation.
- Square Root Semi Kurtosis.
- First Lower Partial Moment (Omega Ratio).
- Second Lower Partial Moment (Sortino Ratio).
- Value at Risk (VaR).
- Conditional Value at Risk (CVaR).
- Entropic Value at Risk (EVaR).
- Relativistic Value at Risk (RLVaR).
- Tail Gini.
- Worst Case Realization (Minimax).

**Drawdown Risk Measures:**

- Average Drawdown for compounded and uncompounded cumulative returns.
- Ulcer Index for compounded and uncompounded cumulative returns.
- Drawdown at Risk (DaR) for compounded and uncompounded cumulative returns.
- Conditional Drawdown at Risk (CDaR) for compounded and uncompounded cumulative returns.
- Entropic Drawdown at Risk (EDaR) for compounded and uncompounded cumulative returns.
- Relativistic Drawdown at Risk (EDaR) for compounded and uncompounded cumulative returns.
- Maximum Drawdown (Calmar Ratio) for compounded and uncompounded cumulative returns.

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
