# Riskfolio-Lib

**Quantitative Strategic Asset Allocation, Easy for Everyone.**

<a href="https://www.kqzyfj.com/click-101360347-15150084?url=https%3A%2F%2Flink.springer.com%2Fbook%2F9783031843037" target="_blank">
<div>
<img src="https://raw.githubusercontent.com/dcajasn/Riskfolio-Lib/refs/heads/master/docs/source/_static/Button.png" height="40" />
</div>
<br>
</a>

<div class="row">
<img src="https://raw.githubusercontent.com/dcajasn/Riskfolio-Lib/master/docs/source/images/MSV_Frontier.png" height="200">
<img src="https://raw.githubusercontent.com/dcajasn/Riskfolio-Lib/master/docs/source/images/Pie_Chart.png" height="200">
</div>

[![](https://img.shields.io/static/v1?label=Sponsor&message=%E2%9D%A4&logo=GitHub&color=%23fe8e86)](https://github.com/sponsors/dcajasn)

<a href='https://ko-fi.com/B0B833SXD' target='_blank'><img height='36' style='border:0px;height:36px;' src='https://cdn.ko-fi.com/cdn/kofi1.png?v=2' border='0' alt='Buy Me a Coffee at ko-fi.com' /></a>

[![GitHub stars](https://img.shields.io/github/stars/dcajasn/Riskfolio-Lib?color=green)](https://github.com/dcajasn/Riskfolio-Lib/stargazers)
[![Downloads](https://static.pepy.tech/badge/Riskfolio-Lib?left_text=Downloads)](https://pepy.tech/project/Riskfolio-Lib)
[![Downloads](https://static.pepy.tech/personalized-badge/riskfolio-lib?period=month&left_color=grey&right_color=orange&left_text=downloads/month)](https://pepy.tech/project/riskfolio-lib)
[![Documentation Status](https://readthedocs.org/projects/riskfolio-lib/badge/?version=latest)](https://riskfolio-lib.readthedocs.io/en/latest/?badge=latest)
[![GitHub license](https://img.shields.io/github/license/dcajasn/Riskfolio-Lib)](https://github.com/dcajasn/Riskfolio-Lib/blob/master/LICENSE.txt)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/dcajasn/Riskfolio-Lib/HEAD)

[![Star History Chart](https://api.star-history.com/svg?repos=dcajasn/Riskfolio-Lib&type=Timeline)](https://star-history.com/#dcajasn/Riskfolio-Lib&Timeline)

## Description

Riskfolio-Lib is a library for making quantitative strategic asset allocation
or portfolio optimization in Python made in Peru &#x1F1F5;&#x1F1EA;. Its objective is to help students, academics and practitioners to build investment portfolios based on mathematically complex models with low effort. It is built on top of
[CVXPY](https://www.cvxpy.org/) and closely integrated
with [Pandas](https://pandas.pydata.org/) data structures.

Some of key functionalities that Riskfolio-Lib offers:

- Mean Risk and Logarithmic Mean Risk (Kelly Criterion) Portfolio Optimization with 4 objective functions:

    - Minimum Risk.
    - Maximum Return.
    - Maximum Utility Function.
    - Maximum Risk Adjusted Return Ratio.

- Mean Risk and Logarithmic Mean Risk (Kelly Criterion) Portfolio Optimization with 24 convex risk measures:

    **Dispersion Risk Measures:**

    - Standard Deviation.
    - Square Root Kurtosis.
    - Mean Absolute Deviation (MAD).
    - Gini Mean Difference (GMD).
    - Conditional Value at Risk Range.
    - Tail Gini Range.
    - Entropic Value at Risk Range.
    - Relativistic Value at Risk Range.
    - Range.
    &nbsp;
    
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
    &nbsp;
    
    **Drawdown Risk Measures:**

    - Average Drawdown for uncompounded cumulative returns.
    - Ulcer Index for uncompounded cumulative returns.
    - Conditional Drawdown at Risk (CDaR) for uncompounded cumulative returns.
    - Entropic Drawdown at Risk (EDaR) for uncompounded cumulative returns.
    - Relativistic Drawdown at Risk (RLDaR) for uncompounded cumulative returns.
    - Maximum Drawdown (Calmar Ratio) for uncompounded cumulative returns.

- Risk Parity Portfolio Optimization with 20 convex risk measures:

    **Dispersion Risk Measures:**

    - Standard Deviation.
    - Square Root Kurtosis.
    - Mean Absolute Deviation (MAD).
    - Gini Mean Difference (GMD).
    - Conditional Value at Risk Range.
    - Tail Gini Range.
    - Entropic Value at Risk Range.
    - Relativistic Value at Risk Range.
    &nbsp;

    **Downside Risk Measures:**

    - Semi Standard Deviation.
    - Square Root Semi Kurtosis.
    - First Lower Partial Moment (Omega Ratio)
    - Second Lower Partial Moment (Sortino Ratio)
    - Conditional Value at Risk (CVaR).
    - Tail Gini.
    - Entropic Value at Risk (EVaR).
    - Relativistic Value at Risk (RLVaR).
    &nbsp;
    
    **Drawdown Risk Measures:**

    - Ulcer Index for uncompounded cumulative returns.
    - Conditional Drawdown at Risk (CDaR) for uncompounded cumulative returns.
    - Entropic Drawdown at Risk (EDaR) for uncompounded cumulative returns.
    - Relativistic Drawdown at Risk (RLDaR) for uncompounded cumulative returns.

- Hierarchical Clustering Portfolio Optimization: Hierarchical Risk Parity (HRP) and Hierarchical Equal Risk Contribution (HERC) with 35 risk measures using naive risk parity:

    **Dispersion Risk Measures:**

    - Standard Deviation.
    - Variance.
    - Square Root Kurtosis.
    - Mean Absolute Deviation (MAD).
    - Gini Mean Difference (GMD).
    - Value at Risk Range.
    - Conditional Value at Risk Range.
    - Tail Gini Range.
    - Entropic Value at Risk Range.
    - Relativistic Value at Risk Range.
    - Range.
    &nbsp;
    
    **Downside Risk Measures:**

    - Semi Standard Deviation.
    - Fourth Root Semi Kurtosis.
    - First Lower Partial Moment (Omega Ratio).
    - Second Lower Partial Moment (Sortino Ratio).
    - Value at Risk (VaR).
    - Conditional Value at Risk (CVaR).
    - Tail Gini.
    - Entropic Value at Risk (EVaR).
    - Relativistic Value at Risk (RLVaR).
    - Worst Case Realization (Minimax).
    &nbsp;
    
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
- Portfolio optimization with constraints based on graph information.
- Portfolio optimization with inequality constraints on risk contributions for variance.
- Portfolio optimization with inequality constraints on factor risk contributions for variance.
- Tools to build efficient frontier for 24 convex risk measures.
- Tools to build linear constraints on assets, asset classes and risk factors.
- Tools to build views on assets and asset classes.
- Tools to build views on risk factors.
- Tools to build risk contribution constraints per asset classes.
- Tools to build risk contribution constraints per risk factor using explicit risk factors and principal components.
- Tools to build bounds constraints for Hierarchical Clustering Portfolios.
- Tools to calculate risk measures.
- Tools to calculate risk contributions per asset.
- Tools to calculate risk contributions per risk factor.
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


## Choosing a Solver

Due to Riskfolio-Lib is based on CVXPY, Riskfolio-Lib can use the same solvers available for CVXPY. The list of solvers compatible with CVXPY is available in [Choosing a solver](https://www.cvxpy.org/tutorial/advanced/index.html#choosing-a-solver) section of CVXPY's documentation. However, to select an adequate solver for each risk measure we can use the following table that specifies which type of programming technique is used to model each risk measure.

| Risk Measure                          | LP | QP | SOCP | SDP | EXP | POW |
|---------------------------------------|----|----|------|-----|-----|-----|
| Variance (MV)                         |    |    | X    | X*  |     |     |
| Mean Absolute Deviation (MAD)         | X  |    |      |     |     |     |
| Gini Mean Difference (GMD)            |    |    |      |     |     | X** |
| Semi Variance (MSV)                   |    |    | X    |     |     |     |
| Kurtosis (KT)                         |    |    |      | X   |     |     |
| Semi Kurtosis (SKT)                   |    |    |      | X   |     |     |
| First Lower Partial Moment (FLPM)     | X  |    |      |     |     |     |
| Second Lower Partial Moment (SLPM)    |    |    | X    |     |     |     |
| Conditional Value at Risk (CVaR)      | X  |    |      |     |     |     |
| Tail Gini (TG)                        |    |    |      |     |     | X** |
| Entropic Value at Risk (EVaR)         |    |    |      |     | X   |     |
| Relativistic Value at Risk (RLVaR)    |    |    |      |     |     | X** |
| Worst Realization (WR)                | X  |    |      |     |     |     |
| CVaR Range (CVRG)                     | X  |    |      |     |     |     |
| Tail Gini Range (TGRG)                |    |    |      |     |     | X** |
| EVaR Range (EVRG)                     |    |    |      |     | X   |     |
| RLVaR Range (RVRG)                    |    |    |      |     |     | X** |
| Range (RG)                            | X  |    |      |     |     |     |
| Average Drawdown (ADD)                | X  |    |      |     |     |     |
| Ulcer Index (UCI)                     |    |    | X    |     |     |     |
| Conditional Drawdown at Risk (CDaR)   | X  |    |      |     |     |     |
| Entropic Drawdown at Risk (EDaR)      |    |    |      |     | X   |     |
| Relativistic Drawdown at Risk (RLDaR) |    |    |      |     |     | X** |
| Maximum Drawdown (MDD)                | X  |    |      |     |     |     |

(*) When SDP graph theory constraints are included. In the case of integer programming graph theory constraints, the model assume the SOCP formulation.

(**) For these models is highly recommended to use MOSEK as solver, due to in some cases CLARABEL cannot find a solution and SCS takes too much time to solve them.

LP - Linear Programming refers to problems with a linear objective function and linear constraints.

QP - Quadratic Programming refers to problems with a quadratic objective function and linear constraints.

SOCP - Second Order Cone Programming refers to problems with second-order cone constraints.

SDP - Semidefinite Programming refers to problems with positive semidefinite constraints.

EXP - refers to problems with exponential cone constraints.

POW - refers to problems with 3-dimensional power cone constraints.


## Dependencies

Riskfolio-Lib supports Python 3.9 or higher.

Installation requires:
- [numpy](http://www.numpy.org/) >= 1.24.0
- [scipy](https://www.scipy.org/) >= 1.10.0
- [pandas](https://pandas.pydata.org/) >= 2.0.0
- [matplotlib](https://matplotlib.org/) >= 3.8.0
- [clarabel](https://oxfordcontrol.github.io/ClarabelDocs/stable/) >= 0.6.0
- [cvxpy](https://www.cvxpy.org/) >= 1.5.2
- [scikit-learn](https://scikit-learn.org/stable/) >= 1.3.0
- [statsmodels](https://www.statsmodels.org/) >= 0.13.5
- [arch](https://bashtage.github.io/arch/) >= 7.0
- [xlsxwriter](https://xlsxwriter.readthedocs.io) >= 3.1.2
- [networkx](https://networkx.org) >= 3.0
- [astropy](https://www.astropy.org) >= 5.1
- [pybind11](https://pybind11.readthedocs.io/en/stable/) >= 2.10.1

## Installation

The latest stable release (and older versions) can be installed from PyPI:

    pip install riskfolio-lib

## Citing

If you use Riskfolio-Lib for published work, please use the following BibTeX entry:

```
@misc{riskfolio,
      author = {Dany Cajas},
      title = {Riskfolio-Lib (7.0.1)},
      year  = {2025},
      url   = {https://github.com/dcajasn/Riskfolio-Lib},
      }
```
 
## Development

Riskfolio-Lib development takes place on Github: https://github.com/dcajasn/Riskfolio-Lib


## Consulting Fees

Riskfolio-Lib is an open-source project, but since it's a project that is not financed for any institution, I started charging for consultancies that are not related to errors in source code. Our fees are as follows:

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

- __[LinkedIn](https://www.linkedin.com/in/dany-cajas/)__
- __[Gmail](dcajasn@gmail.com)__

You can pay using one of the following channels:

- __[Github Sponsorship](https://github.com/sponsors/dcajasn)__

- <a href='https://ko-fi.com/B0B833SXD' target='_blank'><img height='36' style='border:0px;height:36px;' src='https://cdn.ko-fi.com/cdn/kofi1.png?v=2' border='0' alt='Buy Me a Coffee at ko-fi.com' /></a>

## RoadMap

The plan for this module is to add more functions that will be very useful
to asset managers.

- Add more functions based on suggestion of users.
