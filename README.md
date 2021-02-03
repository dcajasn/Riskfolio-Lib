# Riskfolio-Lib

**Quantitative Strategic Asset Allocation, Easy for Everyone.**

<div class="row">
<img src="https://raw.githubusercontent.com/dcajasn/Riskfolio-Lib/master/docs/source/images/MSV_Frontier.png" height="200">
<img src="https://raw.githubusercontent.com/dcajasn/Riskfolio-Lib/master/docs/source/images/Pie_Chart.png" height="200">
</div>

<a href='https://ko-fi.com/B0B833SXD' target='_blank'><img height='36' style='border:0px;height:36px;' src='https://cdn.ko-fi.com/cdn/kofi1.png?v=2' border='0' alt='Buy Me a Coffee at ko-fi.com' /></a>

[![GitHub stars](https://img.shields.io/github/stars/dcajasn/Riskfolio-Lib?color=green)](https://github.com/dcajasn/Riskfolio-Lib/stargazers)
[![Downloads](https://static.pepy.tech/personalized-badge/riskfolio-lib?period=month&units=none&left_color=grey&right_color=orange&left_text=Downloads/Month)](https://pepy.tech/project/riskfolio-lib)
[![Documentation Status](https://readthedocs.org/projects/riskfolio-lib/badge/?version=latest)](https://riskfolio-lib.readthedocs.io/en/latest/?badge=latest)
[![GitHub license](https://img.shields.io/github/license/dcajasn/Riskfolio-Lib)](https://github.com/dcajasn/Riskfolio-Lib/blob/master/LICENSE.txt)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/dcajasn/Riskfolio-Lib/HEAD)


## Description

Riskfolio-Lib is a library for making quantitative strategic asset allocation
or portfolio optimization in Python made in Peru &#x1F1F5;&#x1F1EA;. Its objective is to help students, academics and practitioners to build investment portfolios based on mathematically complex models with low effort. It is built on top of
[cvxpy](https://www.cvxpy.org/) and closely integrated
with [pandas](https://pandas.pydata.org/) data structures.

Some of key functionalities that Riskfolio-Lib offers:

* Mean Risk Portfolio optimization with 4 objective functions:

    * Minimum Risk.
    * Maximum Return.
    * Maximum Utility Function.
    * Maximum Risk Adjusted Return Ratio.

* Mean Risk Portfolio optimization with 13 convex risk measures:

    * Standard Deviation.
    * Semi Standard Deviation.
    * Mean Absolute Deviation (MAD).
    * First Lower Partial Moment (Omega Ratio)
    * Second Lower Partial Moment (Sortino Ratio)
    * Conditional Value at Risk (CVaR).
    * Entropic Value at Risk (EVaR).
    * Worst Case Realization (Minimax Model)
    * Maximum Drawdown (Calmar Ratio)
    * Average Drawdown
    * Conditional Drawdown at Risk (CDaR).
    * Entropic Drawdown at Risk (EDaR).
    * Ulcer Index.

* Risk Parity Portfolio optimization with 10 convex risk measures:

    * Standard Deviation.
    * Semi Standard Deviation.
    * Mean Absolute Deviation (MAD).
    * First Lower Partial Moment (Omega Ratio)
    * Second Lower Partial Moment (Sortino Ratio)
    * Conditional Value at Risk (CVaR).
    * Entropic Value at Risk (EVaR).
    * Conditional Drawdown at Risk (CDaR).
    * Entropic Drawdown at Risk (EDaR).
    * Ulcer Index.

* Worst Case Mean Variance Portfolio optimization.
* Portfolio optimization with Black Litterman model.
* Portfolio optimization with Risk Factors model.
* Portfolio optimization with constraints on tracking error and turnover.
* Portfolio optimization with short positions and leveraged portfolios.
* Tools to build efficient frontier for 13 risk measures.
* Tools to build linear constraints on assets, asset classes and risk factors.
* Tools to build views on assets and asset classes.
* Tools to calculate risk measures.
* Tools to calculate risk contributions per asset.
* Tools to calculate uncertainty sets for mean vector and covariance matrix.
* Tools to estimate loadings matrix (Stepwise Regression and Principal Components Regression).
* Tools to visualizing portfolio properties and risk measures.
* Tools to build reports on Jupyter Notebook and Excel. 
* Option to use commercial optimization solver like MOSEK or GUROBI for large scale problems. 


## Documentation

Online documentation is available at [Documentation](https://riskfolio-lib.readthedocs.io/en/latest/).

The docs include a [tutorial](https://riskfolio-lib.readthedocs.io/en/latest/examples.html)
with examples that shows the capacities of Riskfolio-Lib.


## Dependencies

Riskfolio-Lib supports Python 3.7+.

Installation requires:
* [numpy](http://www.numpy.org/) >= 1.17.0
* [scipy](https://www.scipy.org/) >= 1.1.0
* [pandas](https://pandas.pydata.org/) >= 1.0.0
* [matplotlib](https://matplotlib.org/) >= 3.3.0
* [cvxpy](https://www.cvxpy.org/) >= 1.0.15
* [scikit-learn](https://scikit-learn.org/stable/) >= 0.22.0
* [statsmodels](https://www.statsmodels.org/) >= 0.10.1
* [arch](https://bashtage.github.io/arch/) >= 4.15
* [xlsxwriter](https://xlsxwriter.readthedocs.io) >= 1.3.7


## Installation

The latest stable release (and older versions) can be installed from PyPI:

    pip install riskfolio-lib

 
## Development

Riskfolio-Lib development takes place on Github: https://github.com/dcajasn/Riskfolio-Lib

## RoadMap

The plan for this module is to add more functions that will be very useful
to asset managers.

* Add Black Litterman for factors models.
* Add functions to estimate Duration, Convexity, Key Rate Durations and Convexities of bonds without embedded options (for loadings matrix).
* Add more functions based on suggestion of users.