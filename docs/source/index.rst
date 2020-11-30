#############
Riskfolio-Lib
#############
**Quantitative Strategic Asset Allocation, easy for you.**

.. image:: images/MSV_Frontier.png
    :width: 45%
    
.. image:: images/Pie_Chart.png
    :width: 45%
    

Description
===========

Riskfolio-Lib is a library for making quantitative strategic asset allocation
or portfolio optimization in Python. It is built on top of
`CVXPY <https://www.cvxpy.org/>`_ and closely integrated
with `pandas <https://pandas.pydata.org/>`_ data structures.

Some of key functionality that Riskfolio-Lib offers:

- Portfolio optimization with 4 objective functions (Minimum Risk, Maximum Return, Maximum Risk Adjusted Return Ratio and Maximum Utility Function)
- Portfolio optimization with 10 convex risk measures (Std. Dev., MAD, CVaR, Maximum Drawdown, among others)
- Risk Parity Portfolio optimization with 7 convex risk measures (Std. Dev., MAD, CVaR, Maximum Drawdown, among others)
- Portfolio optimization with Black Litterman model.
- Portfolio optimization with Risk Factors model.
- Portfolio optimization with constraints on tracking error and turnover.
- Portfolio optimization with short positions and leveraged portfolios.
- Tools for build efficient frontier for 10 risk measures.
- Tools for build linear constraints on assets, asset classes and risk factors.
- Tools for build views on assets and asset classes.
- Tools for calculate risk measures.
- Tools for calculate risk contributions per asset.
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

The plan for this module is to add more functions that will be very useful
to asset managers.

* Mean Entropic Risk Optimization Portfolios.
* Mean Risk Worst Case Optimization (Min Max):

    - Box and elipsoidal constraints for mean and covariance matrix.
    - Worst covariance and mean estimation using bootstrapping.
    - Worst covariance and mean estimation using percentage change.
* Add functions to estimate Duration, Convexity, Key Rate Durations and Convexities of bonds without embedded options (for loadings matrix).
* Add more functions based on suggestion of users.
