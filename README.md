# Riskfolio-Lib

**Quantitative Strategic Asset Allocation, easy for you.**

<div class="row">
<img src="https://raw.githubusercontent.com/dcajasn/Riskfolio-Lib/master/docs/source/images/MSV_Frontier.png" height="200">
<img src="https://raw.githubusercontent.com/dcajasn/Riskfolio-Lib/master/docs/source/images/Pie_Chart.png" height="200">
</div>


## Description

Riskfolio-Lib is a library for making quantitative strategic asset allocation
or portfolio optimization in Python. It is built on top of
[cvxpy](https://www.cvxpy.org/) and closely integrated
with [pandas](https://pandas.pydata.org/) data structures.

Some of key functionality that Riskfolio-Lib offers:

* Portfolio optimization with 4 objective functions (Minimum Risk, Maximum Return, Maximum Risk Adjusted Return Ratio and Maximum Utility Function)
* Portfolio optimization with 10 convex risk measures (Std. Dev., MAD, CVaR, Maximum Drawdown, among others)
* Risk Parity Portfolio optimization with 7 convex risk measures (Std. Dev., MAD, CVaR, Maximum Drawdown, among others)
* Portfolio optimization with Black Litterman model.
* Portfolio optimization with Risk Factors model.
* Portfolio optimization with constraints on tracking error and turnover.
* Portfolio optimization with short positions and leveraged portfolios.
* Tools for build efficient frontier for 10 risk measures.
* Tools for build linear constraints on assets, asset classes and risk factors.
* Tools for build views on assets and asset classes.
* Tools for calculate risk measures.
* Tools for calculate risk contributions per asset.
* Tools for estimate loadings matrix (Stepwise Regression and Principal Components Regression).
* Tools for visualizing portfolio properties and risk measures.


## Documentation

Online documentation is available at [Documentation](https://riskfolio-lib.readthedocs.io/en/latest/).

The docs include a [tutorial](https://riskfolio-lib.readthedocs.io/en/latest/examples.html)
with examples that shows the capacities of Riskfolio-Lib.


## Dependencies

Riskfolio-Lib supports Python 3.6+.

Installation requires:
* [numpy](http://www.numpy.org/) >= 1.17.0
* [scipy](https://www.scipy.org/) >= 1.0.1
* [pandas](https://pandas.pydata.org/) >= 1.0.0
* [matplotlib](https://matplotlib.org/) >= 3.0.0
* [cvxpy](https://www.cvxpy.org/) >= 1.0.15
* [scikit-learn](https://scikit-learn.org/stable/) >= 0.22.0
* [statsmodels](https://www.statsmodels.org/) >= 0.10.1

## Installation

The latest stable release (and older versions) can be installed from PyPI:

    pip install riskfolio-lib

 
## Development

Riskfolio-Lib development takes place on Github: https://github.com/dcajasn/Riskfolio-Lib


## Module Plans

The plan for this module is to add more functions that will be very useful
to asset managers.

* Mean Entropic Risk Optimization Portfolios.
* Mean Risk Worst Case Optimization (Min Max):
    - Box and elipsoidal constraints for mean and covariance matrix.
    - Worst covariance and mean estimation using bootstrapping.
    - Worst covariance and mean estimation using percentage change.
* Add functions to estimate Duration, Convexity, Key Rate Durations and Convexities of bonds without embedded options (for loadings matrix).
* Add more functions based on suggestion of users.