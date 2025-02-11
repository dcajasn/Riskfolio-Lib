#########
Changelog
#########

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

Version 7.0.0
=============

- Add two new convex risk measures: EVaR Range and RLVaR Range, to the Portfolio object.
- Add three new risk measures: VaR Range, EVaR Range and RLVaR Range, to the HCPortfolio object.
- Add the generalization of risk parity for variance through inequality constraints on the risk contributions of assets to the Portfolio object.
- Add the generalization of factor risk parity for variance through inequality constraints on the risk contributions of risk factors to the Portfolio object.
- Add a function to calculate the Brinson Performance Attribution per class and aggregate.
- Add a plot function to show the Brinson Performance Attribution.
- Add functions to calculate the VaR Range, EVaR Range and RLVaR Range.
- Update plot functions to consider EVaR Range and RLVaR Range.
- Update duplication, elimination and summation matrices functions to consider or not the diagonal of the symmetric matrix.


Version 6.3.0
=============

- Add new functions to calculate the number of effective assets (NEA) and the average centrality of the portfolio.
- Add the possibility to use neighborhood and cluster network constraints at the same time.
- Fixed some bugs in HRP and HERC when we add constraints.
- Fixed a bug in the duplication_summation_matrix.
- Fixed tight layout in plot functions that uses multiple axes.

Version 6.2.0
=============

- Improvement in calculation speed of duplication_matrix, duplication_elimination_matrix and duplication_summation_matrix functions using a vectorized formula.
- Fixed formulation of risk parity with risk factors model that produced incorrect results when using the MOSEK solver.
- Fixed some bugs in PlotFunctions module.
- Fixed some bugs in HCPortfolio related to custom_mu vector and use of Kurtosis and Semi Kurtosis as risk measures.
- Standardized the way additional parameters to estimate mean vector and covariance matrix are entered.

Version 6.1.0
=============

- Implements standarized silhouette score to determine the optimal number of clusters.
- Fix plot_clusters function to plot clusters and heatmap in same order of codependence matrix. Originally it plots the codependece matrix with axis x inverted.

Version 6.0.0
=============

- Implements risk parity optimization based on explicit risk factors and principal components.
- Implements new formulations of Gini Mean Difference, Tail Gini, Range, CVaR Range and Tail Gini Range that improves speed compared to formulations based on the owa portfolio model.
- Improves the calculation of elliptical uncertainty sets for worst case optimization.
- Add new functions that allow us to calculate the risk contribution per explicit risk factors and principal components.
- Add new functions that allow us to plot the risk contribution per explicit risk factors and principal components.

Version 5.0.0
=============

- Implements new kind of constraints that incorporates the information from networks like the Minimum Spanning Tree and Maximally Filtered Graph into the portfolio optimization models: return-risk portfolio, owa portfolio and worst case portfolio.
- Implements new kind of constraints that incorporates the information from dendrograms  into the portfolio optimization models: return-risk portfolio, owa portfolio and worst case portfolio.
- Improves the speed of several functions using the c++ linear algebra library Eigen and c++ eigenvalues library Spectra.
- Add new functions that allow us to plot the relationship between graphs and asset allocation.
- Add new functions that allow us to create constraints based on graphs information.
- Add a new example about applications of networks and dendrograms constraints in portfolio optimization problems.
- Fixed some errors related to HCPortfolio with constraints.
- Fixed some errors in some plots.

Version 4.4.0
=============

- Implements the approximate Kurtosis model through sum of squared quadratic forms for large scale kurtosis optimization.
- Add the block vectorization operator.

Version 4.3.0
=============

- Implements custom constraints for the Relaxed Risk Parity portfolio model.
- Add three new methods to estimate the mean vector: James-Stein, Bayes-Stein and BOP.

Version 4.2.0
=============

- Implements constraints for the Hierarchical Equal Risk Contribution (HERC) and Nested Clustered Optimization (NCO) portfolio models.
- Add the option to show risk contributions as a percentage of total risk in risk contribution plot.
- Repairs some bugs.

Version 4.1.0
=============

- Implements the Relativistic Value at Risk and Relativistic Drawdown at Risk portfolio models.
- Implements the Higher L-moments portfolio model function as an special case of OWA portfolio.
- Adds functions to calculate L-moments.
- Adds a function to calculate risk contribution constraints on asset classes.
- Repairs some bugs.

Version 4.0.0
=============

- Implements Kurtosis and Semi Kurtosis portfolio models based on parametric approach.
- Implements new c++ based functions to speed up kurtosis model calculations.
- Repairs some bugs.

Version 3.3.0
=============

- Adds Kendall Tau and Gerber statistic as options for codependence matrix in HCPortfolio object.
- Adds Gerber statistic as an option for covariance matrix estimator in Portfolio and HCPortfolio objects.

Version 3.2.0
=============

- Implements reformulations of portfolio models based on drawdowns to speed up calculations.
- Adds some tests for portfolio object and hcportfolio object.

Version 3.1.0
=============

- Implements a reformulation of OWA portfolio optimization to speed up calculations.

Version 3.0.0
=============

- Implements 5 additional risk measures for mean risk model: Gini Mean Difference, Tail Gini, Range, CVaR range and Tail Gini range.
- Implements 4 additional risk measures for risk parity model: Gini Mean Difference, Tail Gini, CVaR range and Tail Gini range.
- Implements the OWA Portfolio Optimization model for custom vector of weights and a module to build OWA weights for some special cases.
- Implements a function to plot range risk measures.
- Adds the option to use Graphical Lasso, j-Logo, denoising and detoning covariance estimates.


Version 2.0.0
=============

- Implement Nested Clustered Optimization (NCO) model with four objective functions.
- Implements the Relaxed Risk Parity model.
- Implements the Risk Budgeting approach for Risk Parity Portfolios with constraints.
- Adds the option to use custom covariance in Hierarchical Clustering Portfolios.

Version 1.0.0
=============

- Redesigns of Riskfolio-Lib interface (Only import riskfolio for all functions).
- Implements Hierarchical Risk Parity (HRP) model with constraints on assets' weights.
- Implements a function that helps to build constraints for the HRP model.
- Implements the Direct Bubble Hierarchical Tree (DBHT) linkage method for HRP and HERC models.
- Implements a function that plots relationship among assets in a network using Minimum Spanning Tree (MST) and Planar Maximally Filtered Graph (PMFG).
- Adds two new codependence measures: mutual information and lower tail dependence index.


Version 0.4.0
=============

- Implements Hierarchical Equal Risk Contribution with equally weights within clusters (HERC2).
- Implements a function that help us to discretize portfolio weights into number of shares given an investment amount.
- Implements the option to select the method to estimate covariance in HRP, HERC and HERC2.
- Adds the option to add constraints on the number of assets and the number of effective assets.
- Fixes an error in two_diff_gap_stat() when number of assets is too small.
- Fixes an error on forward_regression() and backward_regression() when there is no significant feature in regression modes using p-value criterion.
- Adds an example that shows how to build HERC2 portfolios.
- Adds an example that shows how to build constraints on the number of assets and number of effective assets.


Version 0.3.0
=============

- Implements Hierarchical Risk Parity (HRP) and Hierarchical Equal Risk Parity (HERC).
- Implements the function plot_clusters() and plot_dendrogram() that help us to identify clusters based on a distance correlation metric.
- Implements the function assets_clusters() that help us to create asset classes based on hierarchical clusters.
- Adds an example that shows how to build Hierarchical Risk Parity portfolios.
- Adds an example that shows how to build Hierarchical Equal Risk Parity portfolios.


Version 0.2.0
=============

- Implements Logarithmic Mean Risk (Kelly Criterion) Portfolio Optimization models.
- Implements the function plot_bar() that help us to plot portfolios with negative weights.
- Adds the option to build dollar neutral portfolios.
- Adds an example that shows how to build Logarithmic Mean Risk (Kelly Criterion) portfolios.
- Adds an example that shows how to build dollar neutral portfolios.


Version 0.1.5
=============

- Adds the option to add a constraint on minimum portfolio return.
- Adds an example of how to add constraints on portfolio return and risk measures.


Version 0.1.4
=============

- Adds Black Litterman with factors in two flavors: Black Litterman Bayesian model and Augmented Black Litterman model.
- Implements factors_views, a function that allows to design views on risk factors for Black Litterman with factors.
- Repairs some bugs.


Version 0.1.2
=============

- Adds Entropic Drawdown at Risk for Mean Risk Portfolio Optimization and Risk Parity Portfolio Optimization.
- Repairs some bugs.


Version 0.1.1
=============

- Repairs some bugs in Portfolio related to Semi Variance and UCI.
- Implements an option to annualize returns and risk in plot_frontier, Jupyter Notebook and Excel reports.
- Adds examples using Vectorbt for Backtesting and MOSEK for large scale problems.


Version 0.1.0
=============

- Repairs some bugs in RiskFunctions.
- Implements the Reports module that helps to build reports on Jupyter Notebook and Excel.
- Implements plot_table, a function that resume some indicators of a portfolio.
- Adds Entropic Value at Risk for Mean Risk Portfolio Optimization and Risk Parity Portfolio Optimization.


Version 0.0.7
=============

- Implements normal assumption method to estimate box and elliptical uncertainty sets for Worst Case Optimization.
- Implements elliptical uncertainty sets for covariance matrix.
- Adds Ulcer Index for Mean Risk Portfolio Optimization and Risk Parity Portfolio Optimization.
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
- Adds Akaike information criterion, Schwarz information criterion, R squared and adjusted R squared feature selection criterions in stepwise regression.


Version 0.0.3
=============

- Repairs some bugs.
- Implements an option for building constraints common for all assets classes.


Version 0.0.2
=============

- Repairs some bugs.


Version 0.0.1
=============

- Implements robust and ewma estimates.
- Implements Black Litterman model and risk factors models.
- Implements mean risk optimization with 10 risk measures.
