#############
Riskfolio-Lib
#############

.. meta::
   :description lang=en: Portfolio Optimization and Quantitative Strategic Asset Allocation in Python.
   :description lang=es: Optimización de Portafolios y Asignación Estratégica de Activos con Python.
   :keywords lang=en: portfolio optimization python, portfolio optimization, Markowitz portfolio optimization, cvar portfolio optimization, asset allocation, strategic asset allocation
   :keywords lang=es: optimización de portafolios python, optimización de portafolios, optimización de portafolios de Markowitz, optimización de portafolios con cvar, asignación de activos, asignación estratégica de activos


**Quantitative Strategic Asset Allocation, Easy for Everyone.**

.. image:: images/MSV_Frontier.png
    :width: 45%
    
.. image:: images/Pie_Chart.png
    :width: 45%

.. image:: https://img.shields.io/static/v1?label=Sponsor&message=%E2%9D%A4&logo=GitHub&color=%23fe8e86
 :target: https://github.com/sponsors/dcajasn

.. raw:: html
   
   <br>

.. raw:: html

    <a href='https://ko-fi.com/B0B833SXD' target='_blank'><img height='36'style='border:0px;height:36px;' src='https://cdn.ko-fi.com/cdn/kofi1.png?v=2' border='0' alt='Buy Me a Coffee at ko-fi.com' /></a><br>


.. image:: https://img.shields.io/github/stars/dcajasn/Riskfolio-Lib?color=green   :alt: GitHub stars   :target: https://github.com/dcajasn/Riskfolio-Lib/stargazers
.. image:: https://static.pepy.tech/personalized-badge/riskfolio-lib?period=month&units=none&left_color=grey&right_color=orange&left_text=Downloads/Month
 :target: https://pepy.tech/project/riskfolio-lib
.. image:: https://readthedocs.org/projects/riskfolio-lib/badge/?version=latest
.. raw:: html

    <a href="https://github.com/dcajasn/Riskfolio-Lib/blob/master/LICENSE.txt"><img alt="GitHub license" src="https://img.shields.io/github/license/dcajasn/Riskfolio-Lib"></a>

.. image:: https://mybinder.org/badge_logo.svg
 :target: https://mybinder.org/v2/gh/dcajasn/Riskfolio-Lib/HEAD


Description
===========

Riskfolio-Lib is a library for making portfolio optimization and quantitative strategic asset allocation in Python made in Peru |:peru:|. Its objective is to help students, academics and practitioners to build investment portfolios based on mathematically complex models with low effort. It is built on top of
`CVXPY <https://www.cvxpy.org/>`_ and closely integrated
with `pandas <https://pandas.pydata.org/>`_ data structures.

Some of key functionalities that Riskfolio-Lib offers:

- Mean Risk and Logarithmic Mean Risk (Kelly Criterion) Portfolio Optimization with 4 objective functions:

    - Minimum Risk.
    - Maximum Return.
    - Maximum Utility Function.
    - Maximum Risk Adjusted Return Ratio.

- Mean Risk and Logarithmic Mean Risk (Kelly Criterion) Portfolio Optimization with 22 convex risk measures:

    **Dispersion Risk Measures:**

    - Standard Deviation.
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
    - Conditional Value at Risk (CVaR).
    - Tail Gini.
    - Entropic Value at Risk (EVaR).
    - Relativistic Value at Risk (RLVaR).
    - Worst Case Realization (Minimax).

    **Drawdown Risk Measures:**

    - Average Drawdown for uncompounded cumulative returns.
    - Ulcer Index for uncompounded cumulative returns.
    - Conditional Drawdown at Risk (CDaR) for uncompounded cumulative returns.
    - Entropic Drawdown at Risk (EDaR) for uncompounded cumulative returns.
    - Relativistic Drawdown at Risk (RLDaR) for uncompounded cumulative returns.
    - Maximum Drawdown (Calmar Ratio) for uncompounded cumulative returns.

- Risk Parity Portfolio Optimization with 18 convex risk measures:

    **Dispersion Risk Measures:**

    - Standard Deviation.
    - Square Root Kurtosis.
    - Mean Absolute Deviation (MAD).
    - Gini Mean Difference (GMD).
    - Conditional Value at Risk Range.
    - Tail Gini Range.

    **Downside Risk Measures:**

    - Semi Standard Deviation.
    - Square Root Semi Kurtosis.
    - First Lower Partial Moment (Omega Ratio)
    - Second Lower Partial Moment (Sortino Ratio)
    - Conditional Value at Risk (CVaR).
    - Tail Gini.
    - Entropic Value at Risk (EVaR).
    - Relativistic Value at Risk (RLVaR).

    **Drawdown Risk Measures:**

    - Ulcer Index for uncompounded cumulative returns.
    - Conditional Drawdown at Risk (CDaR) for uncompounded cumulative returns.
    - Entropic Drawdown at Risk (EDaR) for uncompounded cumulative returns.
    - Relativistic Drawdown at Risk (RLDaR) for uncompounded cumulative returns.

- Hierarchical Clustering Portfolio Optimization: Hierarchical Risk Parity (HRP) and Hierarchical Equal Risk Contribution (HERC) with 32 risk measures using naive risk parity:

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
    - Relativistic Drawdown at Risk (RLDaR) for compounded and uncompounded cumulative returns.
    - Maximum Drawdown (Calmar Ratio) for compounded and uncompounded cumulative returns.

- Nested Clustered Optimization (NCO) with four objective functions and the available risk measures to each objective:

    - Minimum Risk.
    - Maximum Return.
    - Maximum Utility Function.
    - Equal Risk Contribution.

- Worst Case Mean Variance Portfolio Optimization.
- Relaxed Risk Parity Portfolio Optimization.
- Ordered Weighted Averaging (OWA) Portfolio Optimization.
- Portfolio optimization with Black Litterman model.
- Portfolio optimization with Risk Factors model.
- Portfolio optimization with Black Litterman Bayesian model.
- Portfolio optimization with Augmented Black Litterman model.
- Portfolio optimization with constraints on tracking error and turnover.
- Portfolio optimization with short positions and leveraged portfolios.
- Portfolio optimization with constraints on number of assets and number of effective assets.
- Tools to build efficient frontier for 22 convex risk measures.
- Tools to build linear constraints on assets, asset classes and risk factors.
- Tools to build views on assets and asset classes.
- Tools to build views on risk factors.
- Tools to build risk contribution constraints per asset classes.
- Tools to build bounds constraints for Hierarchical Clustering Portfolios.
- Tools to calculate risk measures.
- Tools to calculate risk contributions per asset.
- Tools to calculate uncertainty sets for mean vector and covariance matrix.
- Tools to calculate assets clusters based on codependence metrics.
- Tools to estimate loadings matrix (Stepwise Regression and Principal Components Regression).
- Tools to visualizing portfolio properties and risk measures.
- Tools to build reports on Jupyter Notebook and Excel. 
- Option to use commercial optimization solver like MOSEK or GUROBI for large scale problems.


Consulting Fees
===============

Riskfolio-Lib is an open-source project, however due it's a project that is not financed for any institution, I started charging for consultancies that are not related to errors in source code. Our fees are as follows:

- $ 25 USD (United States Dollars) per question that doesn't require to check code.

- $ 50 USD to check a small size script or code (less than 200 lines of code). The fee of the solution depends on the complexity of the solution:
    - $ 50 USD for simple errors in scripts (modify less than 10 lines of code).
    - For most complex errors the fee depends on the complexity of the solution but the fee is $ 150 USD per hour.

- $ 100 USD to check a medium size script or code (between 201 and 600 lines of code). The fee of the solution depends on the complexity of the solution:
    - $ 50 USD for simple errors in scripts (modify less than 10 lines of code).
    - For most complex errors the fee depends on the complexity of the solution but the fee is $ 150 USD per hour.

- For large size script or code (more than 600 lines of code) the fee is variable depending on the size of the code. The fee of the solution depends on the complexity of the solution:
    - $ 50 USD for simple errors in scripts (modify less than 10 lines of code).
    - For most complex errors the fee depends on the complexity of the solution but the fee is $ 150 USD per hour.

**All consulting must be paid in advance**.

You can contact me through:

- __`LinkedIn <https://www.linkedin.com/in/dany-cajas/>´__
- __`Gmail <dcajasn@gmail.com>`__

You can pay using one of the following channels:

.. image:: https://img.shields.io/static/v1?label=Sponsor&message=%E2%9D%A4&logo=GitHub&color=%23fe8e86
 :target: https://github.com/sponsors/dcajasn
 
.. raw:: html

    <a href='https://ko-fi.com/B0B833SXD' target='_blank'><img height='36'style='border:0px;height:36px;' src='https://cdn.ko-fi.com/cdn/kofi1.png?v=2' border='0' alt='Buy Me a Coffee at ko-fi.com' /></a><br>
    
    
Citing
======

If you use Riskfolio-Lib for published work, please use the following BibTeX entry:

::

    @misc{riskfolio,
          author = {Dany Cajas},
          title = {Riskfolio-Lib (4.4.1)},
          year  = {2023},
          url   = {https://github.com/dcajasn/Riskfolio-Lib},
          }

Contents
========

..  toctree::
    :maxdepth: 1

    Install <install>
    Portfolio Models <portfolio>
    Hierarchical Clustering Models <hcportfolio>
    Parameters Estimation <parameters>
    Constraints Functions <constraints>
    Risk Functions <risk>
    Plot Functions <plot>
    Reports <reports>
    Auxiliary Functions <auxiliary>
    Examples <examples>
    Contributing <contributing>
    Authors <authors>
    License <license>
    Changelog <changelog>

    
Indices and tables
==================

- :ref:`genindex`
- :ref:`modindex`
- :ref:`search`


Module Plans
==================

The plan for this library is to add more functions that will be very useful
for students, academics and practitioners.

- Add more functions based on suggestion of users.
