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

[![Star History Chart](https://api.star-history.com/svg?repos=dcajasn/Riskfolio-Lib&type=Timeline)](https://star-history.com/#dcajasn/Riskfolio-Lib&Timeline)

## Description

Riskfolio-Lib is a library for making quantitative strategic asset allocation
or portfolio optimization in Python made in Peru &#x1F1F5;&#x1F1EA;. Its objective is to help students, academics and practitioners to build investment portfolios based on mathematically complex models with low effort. It is built on top of
[cvxpy](https://www.cvxpy.org/) and closely integrated
with [pandas](https://pandas.pydata.org/) data structures.

Some of key functionalities that Riskfolio-Lib offers:

- Mean Risk and Logarithmic Mean Risk (Kelly Criterion) Portfolio Optimization with 4 objective functions:

    - Minimum Risk.
    - Maximum Return.
    - Maximum Utility Function.
    - Maximum Risk Adjusted Return Ratio.

- Mean Risk and Logarithmic Mean Risk (Kelly Criterion) Portfolio Optimization with 18 convex risk measures:

    **Dispersion Risk Measures:**

    - Standard Deviation.
    - Mean Absolute Deviation (MAD).
    - Gini Mean Difference (GMD).
    - Range.
    - Conditional Value at Risk Range.
    - Tail Gini Range.
    &nbsp;
    
    **Downside Risk Measures:**

    - Semi Standard Deviation.
    - First Lower Partial Moment (Omega Ratio).
    - Second Lower Partial Moment (Sortino Ratio).
    - Conditional Value at Risk (CVaR).
    - Tail Gini.
    - Entropic Value at Risk (EVaR).
    - Worst Case Realization (Minimax).
    &nbsp;
    
    **Drawdown Risk Measures:**

    - Maximum Drawdown (Calmar Ratio) for uncompounded cumulative returns.
    - Average Drawdown for uncompounded cumulative returns.
    - Conditional Drawdown at Risk (CDaR) for uncompounded cumulative returns.
    - Entropic Drawdown at Risk (EDaR) for uncompounded cumulative returns.
    - Ulcer Index for uncompounded cumulative returns.

- Risk Parity Portfolio Optimization with 14 convex risk measures:

    **Dispersion Risk Measures:**

    - Standard Deviation.
    - Mean Absolute Deviation (MAD).
    - Gini Mean Difference (GMD).
    - Conditional Value at Risk Range.
    - Tail Gini Range.
    &nbsp;
    
    **Downside Risk Measures:**

    - Semi Standard Deviation.
    - First Lower Partial Moment (Omega Ratio)
    - Second Lower Partial Moment (Sortino Ratio)
    - Conditional Value at Risk (CVaR).
    - Tail Gini.
    - Entropic Value at Risk (EVaR).
    &nbsp;
    
    **Drawdown Risk Measures:**

    - Conditional Drawdown at Risk (CDaR) for uncompounded cumulative returns.
    - Entropic Drawdown at Risk (EDaR) for uncompounded cumulative returns.
    - Ulcer Index for uncompounded cumulative returns.

- Hierarchical Clustering Portfolio Optimization: Hierarchical Risk Parity (HRP) and Hierarchical Equal Risk Contribution (HERC) with 22 risk measures using naive risk parity:

    **Dispersion Risk Measures:**

    - Standard Deviation.
    - Variance.
    - Mean Absolute Deviation (MAD).
    - Range.
    - Conditional Value at Risk Range.
    - Tail Gini Range.
    &nbsp;
    
    **Downside Risk Measures:**

    - Semi Standard Deviation.
    - First Lower Partial Moment (Omega Ratio).
    - Second Lower Partial Moment (Sortino Ratio).
    - Value at Risk (VaR).
    - Conditional Value at Risk (CVaR).
    - Entropic Value at Risk (EVaR).
    - Tail Gini.
    - Worst Case Realization (Minimax).
    &nbsp;
    
    **Drawdown Risk Measures:**

    - Maximum Drawdown (Calmar Ratio) for compounded and uncompounded cumulative returns.
    - Average Drawdown for compounded and uncompounded cumulative returns.
    - Drawdown at Risk (DaR) for compounded and uncompounded cumulative returns.
    - Conditional Drawdown at Risk (CDaR) for compounded and uncompounded cumulative returns.
    - Entropic Drawdown at Risk (EDaR) for compounded and uncompounded cumulative returns.
    - Ulcer Index for compounded and uncompounded cumulative returns.

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
- Tools to build efficient frontier for 18 risk measures.
- Tools to build linear constraints on assets, asset classes and risk factors.
- Tools to build views on assets and asset classes.
- Tools to build views on risk factors.
- Tools to calculate risk measures.
- Tools to calculate risk contributions per asset.
- Tools to calculate uncertainty sets for mean vector and covariance matrix.
- Tools to calculate assets clusters based on codependence metrics.
- Tools to estimate loadings matrix (Stepwise Regression and Principal Components Regression).
- Tools to visualizing portfolio properties and risk measures.
- Tools to build reports on Jupyter Notebook and Excel. 
- Option to use commercial optimization solver like MOSEK or GUROBI for large scale problems.


## Documentation

Online documentation is available at [Documentation](https://riskfolio-lib.readthedocs.io/en/latest/).

The docs include a [tutorial](https://riskfolio-lib.readthedocs.io/en/latest/examples.html)
with examples that shows the capacities of Riskfolio-Lib.


## Dependencies

Riskfolio-Lib supports Python 3.7+.

Installation requires:
- [numpy](http://www.numpy.org/) >= 1.17.0
- [scipy](https://www.scipy.org/) >= 1.1.0
- [pandas](https://pandas.pydata.org/) >= 1.0.0
- [matplotlib](https://matplotlib.org/) >= 3.3.0
- [cvxpy](https://www.cvxpy.org/) >= 1.0.15
- [scikit-learn](https://scikit-learn.org/stable/) >= 1.0.0
- [statsmodels](https://www.statsmodels.org/) >= 0.10.1
- [arch](https://bashtage.github.io/arch/) >= 4.15
- [xlsxwriter](https://xlsxwriter.readthedocs.io) >= 1.3.7
- [networkx](https://networkx.org) >= 2.5.1
- [astropy](https://www.astropy.org) >= 4.3.1

## Installation

The latest stable release (and older versions) can be installed from PyPI:

    pip install riskfolio-lib

## Citing

If you use Riskfolio-Lib for published work, please use the following BibTeX entrie:

```
@misc{riskfolio,
      author = {Dany Cajas},
      title = {Riskfolio-Lib (3.0.1)},
      year  = {2021},
      url   = {https://github.com/dcajasn/Riskfolio-Lib},
      }
```
 
## Development

Riskfolio-Lib development takes place on Github: https://github.com/dcajasn/Riskfolio-Lib

## RoadMap

The plan for this module is to add more functions that will be very useful
to asset managers.

- Add more functions based on suggestion of users.