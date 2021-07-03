#########
Changelog
#########

.. raw:: html

    <a href='https://ko-fi.com/B0B833SXD' target='_blank'><img height='36'style='border:0px;height:36px;' src='https://cdn.ko-fi.com/cdn/kofi1.png?v=2' border='0' alt='Buy Me a Coffee at ko-fi.com' /></a>

Version 0.4.0
=============

- Implements Hierarchical Equal Risk Parity with equally weights within clusters (HERC2).
- Implements a function that help us to discretize portfolio weights into number of shares given an investment amount.
- Implements the option to select the method to estimate covariance in HRP, HERC and HERC2.
- Add the option to add constraints on the number of assets and the number of effective assets.
- Fix an error in two_diff_gap_stat() when number of assets is too small.
- Fix an error on forward_regression() and backward_regression() when there is no significant feature in regression modes using p-value criterion.
- Add an example that shows how to build HERC2 portfolios.
- Add an example that shows how to build constraints on the number of assets and number of effective assets.


Version 0.3.0
=============

- Implements Hierarchical Risk Parity (HRP) and Hierarchical Equal Risk Parity (HERC).
- Implements the function plot_clusters() and plot_dendrogram() that help us to identify clusters based on a distance correlation metric.
- Implements the function assets_clusters() that help us to create asset classes based on hierarchical clusters.
- Add an example that shows how to build Hierarchical Risk Parity portfolios.
- Add an example that shows how to build Hierarchical Equal Risk Parity portfolios.


Version 0.2.0
=============

- Implements Logarithmic Mean Risk (Kelly Criterion) Portfolio Optimization models.
- Implements the function plot_bar() that help us to plot portfolios with negative weights.
- Add the option to build dollar neutral portfolios.
- Add an example that shows how to build Logarithmic Mean Risk (Kelly Criterion) portfolios.
- Add an example that shows how to build dollar neutral portfolios.


Version 0.1.5
=============

- Add the option to add a constraint on minimum portfolio return.
- Add an example of how to add constraints on portfolio return and risk measures.


Version 0.1.4
=============

- Add Black Litterman with factors in two flavors: Black Litterman Bayesian model and Augmented Black Litterman model.
- Implements factors_views, a function that allows to design views on risk factors for Black Litterman with factors.
- Repairs some bugs.


Version 0.1.2
=============

- Add Entropic Drawdown at Risk for Mean Risk Portfolio Optimization and Risk Parity Portfolio Optimization.
- Repairs some bugs.


Version 0.1.1
=============

- Repairs some bugs in Portfolio related to Semi Variance and UCI.
- Implements an option to annualize returns and risk in plot_frontier, Jupyter Notebook and Excel reports.
- Add examples using Vectorbt for Backtesting and MOSEK for large scale problems.


Version 0.1.0
=============

- Repairs some bugs in RiskFunctions.
- Implements the Reports module that helps to build reports on Jupyter Notebook and Excel.
- Implements plot_table, a function that resume some indicators of a portfolio.
- Add Entropic Value at Risk for Mean Risk Portfolio Optimization and Risk Parity Portfolio Optimization.


Version 0.0.7
=============

- Implements normal assumption method to estimate box and elliptical uncertainty sets for Worst Case Optimization.
- Implements elliptical uncertainty sets for covariance matrix.
- Add Ulcer Index for Mean Risk Portfolio Optimization and Risk Parity Portfolio Optimization.
- Implements functions to calculate Ulcer Index.


Version 0.0.6
=============

- Repairs some bugs.
- Implements bootstrapping methods to estimate box and elliptical uncertainty sets for Worst Case Optimization.
- Implements Worst Case Mean Variance Portfolio Optimization using box and elliptical uncertainty sets.


Version 0.0.5
=============

- Repairs some bugs.
- Implements Risk Parity Portfolio Optimization for 7 convex risk measures.


Version 0.0.4
=============

- Repairs some bugs.
- Update to make it compatible with cvxpy >=1.1.0
- Implements Principal Component Regression for loadings matrix estimation.
- Add Akaike information criterion, Schwarz information criterion, R squared and adjusted R squared feature selection criterions in stepwise regression.


Version 0.0.3
=============

- Repairs some bugs.
- Implements an option for building constraints common for all assets classes.


Version 0.0.2
=============

- Repairs some bugs.


Version 0.0.1
=============

- Implements robust estimates and ewma estimates.
- Implements Black Litterman model and risk factors models.
- Implements mean risk optimization with 10 risk measures.
