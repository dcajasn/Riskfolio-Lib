=============
Riskfolio-Lib
=============

**Quantitative Strategic Asset Allocation, easy for you.**

.. image:: docs/source/images/MSV_Frontier.png
    :width: 45%
.. image:: docs/source/images/Pie_Chart.png
    :width: 45%


Description
-----------

Riskfolio-Lib is a library for making quantitative strategic asset allocation
or portfolio optimization in Python. It is built on top of
`CVXPY <https://www.cvxpy.org/>`_ and closely integrated
with `pandas <https://pandas.pydata.org/>`_ data structures.

Some of key functionality that Riskfolio-Lib offers:

- Portfolio optimization with 4 objective functions (Minimum Risk, Maximum Return, Maximum Risk Adjusted Return Ratio and Maximum Utility Function)
- Portfolio optimization with 10 convex risk measures (Std. Dev., MAD, CVaR, Maximum Drawdown, among others)
- Portfolio optimization with Black Litterman model.
- Portfolio optimization with Risk Factors model.
- Portfolio optimization with constraints on tracking error and turnover.
- Portfolio optimization with short positions and leverage.
- Tools for construct efficient frontier for 10 risk measures.
- Tools for construct linear constraints on assets, asset classes and risk factors.
- Tools for construct views on assets and asset classes.
- Tools for calculate risk measures.
- Tools for visualizing portfolio properties and risk measures.


Documentation
-------------

Online documentation is available at `Documentation https://www.google.com>`_.

The docs include a `tutorial https://www.google.com>`_ with
examples that shows the capacities of Riskfolio-Lib.


Dependencies
------------

Riskfolio-Lib supports Python 3.6+.

Installation requires:

- `numpy <http://www.numpy.org/>`_ >= 1.17.0
- `scipy <https://www.scipy.org/>`_ >= 1.0.1
- `pandas <https://pandas.pydata.org/>`_ >= 1.0.0
- `matplotlib https://matplotlib.org/>`_ >= 3.0.0
- `cvxpy <https://www.cvxpy.org/>`_ >= 1.0.15
- `scikit-learn https://scikit-learn.org/stable/>`_ >= 0.22.0
- `statsmodels (https://www.statsmodels.org/>`_ >= 0.10.1

Installation
------------

The latest stable release (and older versions) can be installed from PyPI:

    pip install riskfolio-lib

 
Development
-----------

Riskfolio-Lib development takes place on Github: https://www.google.com
