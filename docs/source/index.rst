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

[![ko-fi](https://www.ko-fi.com/img/githubbutton_sm.svg)](https://ko-fi.com/B0B833SXD)

[![PyPi](https://img.shields.io/pypi/v/Riskfolio-Lib.svg)]((https://pypi.org/project/Riskfolio-Lib/))
[![Python](https://img.shields.io/pypi/pyversions/Riskfolio-Lib.svg)]((https://pypi.org/project/Riskfolio-Lib/))
[![Documentation Status](https://readthedocs.org/projects/riskfolio-lib/badge/?version=latest)](https://riskfolio-lib.readthedocs.io/en/latest/)
[![GitHub license](https://img.shields.io/github/license/dcajasn/Riskfolio-Lib)](https://github.com/dcajasn/Riskfolio-Lib/blob/master/LICENSE.txt)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/dcajasn/Riskfolio-Lib/HEAD)
[![GitHub stars](https://img.shields.io/github/stars/dcajasn/Riskfolio-Lib?style=social)](https://github.com/dcajasn/Riskfolio-Lib/stargazers)


Description
===========

Riskfolio-Lib is a library for making quantitative strategic asset allocation
or portfolio optimization in Python made in Peru |:peru:|. It is built on top of
`CVXPY <https://www.cvxpy.org/>`_ and closely integrated
with `pandas <https://pandas.pydata.org/>`_ data structures.

Some of key functionalities that Riskfolio-Lib offers:

- Mean Risk Portfolio optimization with 4 objective functions:

    - Minimum Risk.
    - Maximum Return.
    - Maximum Utility Function.
    - Maximum Risk Adjusted Return Ratio.

- Mean Risk Portfolio optimization with 11 convex risk measures:

    - Standard Deviation.
    - Semi Standard Deviation.
    - Mean Absolute Deviation (MAD).
    - First Lower Partial Moment (Omega Ratio)
    - Second Lower Partial Moment (Sortino Ratio)
    - Conditional Value at Risk (CVaR).
    - Worst Case Realization (Minimax Model)
    - Maximum Drawdown (Calmar Ratio).
    - Average Drawdown
    - Conditional Drawdown at Risk (CDaR).
    - Ulcer Index.

- Risk Parity Portfolio optimization with 8 convex risk measures:

    - Standard Deviation.
    - Semi Standard Deviation.
    - Mean Absolute Deviation (MAD).
    - First Lower Partial Moment (Omega Ratio)
    - Second Lower Partial Moment (Sortino Ratio)
    - Conditional Value at Risk (CVaR).
    - Conditional Drawdown at Risk (CDaR).
    - Ulcer Index.

- Worst Case Mean Variance Portfolio optimization.
- Portfolio optimization with Black Litterman model.
- Portfolio optimization with Risk Factors model.
- Portfolio optimization with constraints on tracking error and turnover.
- Portfolio optimization with short positions and leveraged portfolios.
- Tools for build efficient frontier for 11 risk measures.
- Tools for build linear constraints on assets, asset classes and risk factors.
- Tools for build views on assets and asset classes.
- Tools for calculate risk measures.
- Tools for calculate risk contributions per asset.
- Tools for calculate uncertainty sets for mean vector and covariance matrix.
- Tools for estimate loadings matrix (Stepwise Regression and Principal Components Regression).
- Tools for visualizing portfolio properties and risk measures.


Contents
========

..  toctree::
    :maxdepth: 1

    Install <install>
    Portfolio Models <portfolio>
    Parameters Estimation <parameters>
    Constraints Functions <constraints>
    Risk Functions <risk>
    Plot Functions <plot>
    Auxiliary Functions <auxiliary>
    Examples <examples>
    Contributing <contributing>
    Authors <authors>
    License <license>
    Changelog <changelog>

    
Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


Module Plans
==================

The plan for this modules is to add more functions that will be very useful
to asset managers.

* Mean Entropic Risk Optimization Portfolios.
* Add functions to estimate Duration, Convexity, Key Rate Durations and Convexities of bonds without embedded options (for loadings matrix).
* Add more functions based on suggestion of users.
