""""""  #

"""
Copyright (c) 2020-2024, Dany Cajas
All rights reserved.
This work is licensed under BSD 3-Clause "New" or "Revised" License.
License available at https://github.com/dcajasn/Riskfolio-Lib/blob/master/LICENSE.txt
"""

import numpy as np
import pandas as pd
import statsmodels.api as sm
import scipy.stats as st
import sklearn.covariance as skcov
import arch.bootstrap as bs

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from numpy.linalg import inv
from itertools import product

import riskfolio.src.AuxFunctions as af
import riskfolio.src.DBHT as db
import riskfolio.src.GerberStatistic as gs
import riskfolio.external.cppfunctions as cf


__all__ = [
    "mean_vector",
    "covar_matrix",
    "cokurt_matrix",
    "forward_regression",
    "backward_regression",
    "PCR",
    "loadings_matrix",
    "risk_factors",
    "black_litterman",
    "augmented_black_litterman",
    "black_litterman_bayesian",
    "bootstrapping",
    "normal_simulation",
]


def mean_vector(X, method="hist", d=0.94, target="b1"):
    r"""
    Calculate the expected returns vector using the selected method.

    Parameters
    ----------
    X : DataFrame of shape (n_samples, n_assets)
        Assets returns DataFrame, where n_samples is the number of
        observations and n_assets is the number of assets.
    method : str, optional
        The method used to estimate the expected returns.
        The default value is 'hist'. Possible values are:

        - 'hist': use historical estimator.
        - 'ewma1': use ewma with adjust=True. For more information see `EWM <https://pandas.pydata.org/pandas-docs/stable/user_guide/window.html#exponentially-weighted-window>`__.
        - 'ewma2': use ewma with adjust=False. For more information see `EWM <https://pandas.pydata.org/pandas-docs/stable/user_guide/window.html#exponentially-weighted-window>`__.
        - 'JS': James-Stein estimator. For more information see :cite:`b-Meucci2005` and :cite:`b-Feng2016`.
        - 'BS': Bayes-Stein estimator. For more information see :cite:`b-Jorion1986`.
        - 'BOP': BOP estimator. For more information see :cite:`b-Bodnar2019`.

    d : scalar
        The smoothing factor of ewma methods.
        The default is 0.94.

    target : str, optional
        The target mean vector. The default value is 'b1'.
        Possible values are:

        - 'b1': grand mean.
        - 'b2': volatility weighted grand mean.
        - 'b3': mean square error of sample mean.

    Returns
    -------
    mu : 1d-array
        The estimation of expected returns.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    """

    if not isinstance(X, pd.DataFrame):
        raise ValueError("X must be a DataFrame")

    assets = X.columns.tolist()

    if method == "hist":
        mu = np.array(X.mean(), ndmin=2)
    elif method == "ewma1":
        mu = np.array(X.ewm(alpha=1 - d).mean().iloc[-1, :], ndmin=2)
    elif method == "ewma2":
        mu = np.array(X.ewm(alpha=1 - d, adjust=False).mean().iloc[-1, :], ndmin=2)
    elif method == "ewma2":
        mu = np.array(X.ewm(alpha=1 - d, adjust=False).mean().iloc[-1, :], ndmin=2)
    elif method in ["JS", "BS", "BOP"]:
        T, n = np.array(X, ndmin=2).shape
        ones = np.ones((n, 1))
        mu = np.array(X.mean(), ndmin=2).reshape(-1, 1)
        Sigma = np.cov(X, rowvar=False)
        Sigma_inv = np.linalg.inv(Sigma)
        eigvals = np.linalg.eigvals(Sigma)

        # Calculate target vector
        if target == "b1":
            b = ones.T @ mu / n * ones
        elif target == "b2":
            b = ones.T @ Sigma_inv @ mu / (ones.T @ Sigma_inv @ ones) * ones
        elif target == "b3":
            b = np.trace(Sigma) / T * ones

        # Calculate Estimators
        if method == "JS":
            alpha_1 = (
                1
                / T
                * (n * np.mean(eigvals) - 2 * np.max(eigvals))
                / ((mu - b).T @ (mu - b))
            )
            mu = (1 - alpha_1) * mu + alpha_1 * b
        elif method == "BS":
            alpha_1 = (n + 2) / ((n + 2) + T * (mu - b).T @ Sigma_inv @ (mu - b))
            mu = (1 - alpha_1) * mu + alpha_1 * b
        elif method == "BOP":
            alpha_1 = (mu.T @ Sigma_inv @ mu - n / (T - n)) * b.T @ Sigma_inv @ b - (
                mu.T @ Sigma_inv @ b
            ) ** 2
            alpha_1 /= (mu.T @ Sigma_inv @ mu) * (b.T @ Sigma_inv @ b) - (
                mu.T @ Sigma_inv @ b
            ) ** 2
            beta_1 = (1 - alpha_1) * (mu.T @ Sigma_inv @ b) / (mu.T @ Sigma_inv @ mu)
            mu = alpha_1 * mu + beta_1 * b
        mu = mu.T

    mu = pd.DataFrame(np.array(mu, ndmin=2), columns=assets)

    return mu


def covar_matrix(
    X,
    method="hist",
    d=0.94,
    alpha=0.1,
    bWidth=0.01,
    detone=False,
    mkt_comp=1,
    threshold=0.5,
):
    r"""
    Calculate the covariance matrix using the selected method.

    Parameters
    ----------
    X : DataFrame of shape (n_samples, n_assets)
        Assets returns DataFrame, where n_samples is the number of
        observations and n_assets is the number of assets.
    method : str, optional
        The method used to estimate the covariance matrix:
        The default is 'hist'. Possible values are:

        - 'hist': use historical estimates.
        - 'semi': use semi lower covariance matrix.
        - 'ewma1': use ewma with adjust=True. For more information see `EWM <https://pandas.pydata.org/pandas-docs/stable/user_guide/window.html#exponentially-weighted-window>`__.
        - 'ewma2': use ewma with adjust=False. For more information see `EWM <https://pandas.pydata.org/pandas-docs/stable/user_guide/window.html#exponentially-weighted-window>`__.
        - 'ledoit': use the Ledoit and Wolf Shrinkage method.
        - 'oas': use the Oracle Approximation Shrinkage method.
        - 'shrunk': use the basic Shrunk Covariance method.
        - 'gl': use the basic Graphical Lasso Covariance method.
        - 'jlogo': use the j-LoGo Covariance method. For more information see: :cite:`b-jLogo`.
        - 'fixed': denoise using fixed method. For more information see chapter 2 of :cite:`b-MLforAM`.
        - 'spectral': denoise using spectral method. For more information see chapter 2 of :cite:`b-MLforAM`.
        - 'shrink': denoise using shrink method. For more information see chapter 2 of :cite:`b-MLforAM`.
        - 'gerber1': use the Gerber statistic 1. For more information see: :cite:`b-Gerber2021`.
        - 'gerber2': use the Gerber statistic 2. For more information see: :cite:`b-Gerber2021`.

    d : scalar
        The smoothing factor of ewma methods. The default is 0.94.
    alpha : scalar
        The shrfactor of shrunk and shrink method. The default is 0.1.
    bWidth : float
        The bandwidth of the kernel for 'fixed', 'spectral' and 'shrink' methods.
    detone : bool, optional
        If remove the first mkt_comp of correlation matrix for 'fixed', 'spectral'
        and 'shrink' methods. The detone correlation matrix is singular, so it
        cannot be inverted.
    mkt_comp : int, optional
        Number of first components that will be removed using the detone method.
    threshold : float
        Threshold for 'gerber1' and 'gerber2' methods is between 0 and 1.

    Returns
    -------
    cov : nd-array
        The estimation of covariance matrix.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    """

    if not isinstance(X, pd.DataFrame):
        raise ValueError("X must be a DataFrame")

    assets = X.columns.tolist()

    if method == "hist":
        cov = np.cov(X, rowvar=False)
    elif method == "semi":
        T, N = X.shape
        mu = X.mean().to_numpy().reshape(1, -1)
        a = X - np.repeat(mu, T, axis=0)
        a = np.minimum(a, np.zeros_like(a))
        cov = 1 / (T - 1) * a.T @ a
    elif method == "ewma1":
        cov = X.ewm(alpha=1 - d).cov()
        item = cov.iloc[-1, :].name[0]
        cov = cov.loc[(item, slice(None)), :]
    elif method == "ewma2":
        cov = X.ewm(alpha=1 - d, adjust=False).cov()
        item = cov.iloc[-1, :].name[0]
        cov = cov.loc[(item, slice(None)), :]
    elif method == "ledoit":
        lw = skcov.LedoitWolf()
        lw.fit(X)
        cov = lw.covariance_
    elif method == "oas":
        oas = skcov.OAS()
        oas.fit(X)
        cov = oas.covariance_
    elif method == "shrunk":
        sc = skcov.ShrunkCovariance(shrinkage=alpha)
        sc.fit(X)
        cov = sc.covariance_
    elif method == "gl":
        gl = skcov.GraphicalLassoCV()
        gl.fit(X)
        cov = gl.covariance_
    elif method == "jlogo":
        S = np.cov(X, rowvar=False)
        R = np.corrcoef(X, rowvar=False)
        D = np.sqrt(np.clip((1 - R) / 2, a_min=0.0, a_max=1.0))
        (_, _, separators, cliques, _) = db.PMFG_T2s(1 - D**2, nargout=4)
        cov = db.j_LoGo(S, separators, cliques)
        cov = np.linalg.inv(cov)
    elif method in ["fixed", "spectral", "shrink"]:
        cov = np.cov(X, rowvar=False)
        T, N = X.shape
        q = T / N
        cov = af.denoiseCov(
            cov,
            q,
            kind=method,
            bWidth=bWidth,
            detone=detone,
            mkt_comp=int(mkt_comp),
            alpha=alpha,
        )
    elif method == "gerber1":
        cov = gs.gerber_cov_stat1(X, threshold=threshold)
    elif method == "gerber2":
        cov = gs.gerber_cov_stat2(X, threshold=threshold)

    cov = pd.DataFrame(np.array(cov, ndmin=2), columns=assets, index=assets)

    return cov


def cokurt_matrix(
    X,
    method="hist",
    alpha=0.1,
    bWidth=0.01,
    detone=False,
    mkt_comp=1,
):
    r"""
    Calculate the cokurtosis square matrix using the selected method.

    Parameters
    ----------
    X : DataFrame of shape (n_samples, n_assets)
        Assets returns DataFrame, where n_samples is the number of
        observations and n_assets is the number of assets.
    method : str, optional
        The method used to estimate the cokurtosis square matrix:
        The default is 'hist'. Possible values are:

        - 'hist': use historical estimates.
        - 'semi': use semi lower cokurtosis square matrix.
        - 'fixed': denoise using fixed method. For more information see chapter 2 of :cite:`b-MLforAM`.
        - 'spectral': denoise using spectral method. For more information see chapter 2 of :cite:`b-MLforAM`.
        - 'shrink': denoise using shrink method. For more information see chapter 2 of :cite:`b-MLforAM`.
    bWidth : float
        The bandwidth of the kernel for 'fixed', 'spectral' and 'shrink' methods.
    detone : bool, optional
        If remove the first mkt_comp of correlation matrix for 'fixed', 'spectral'
        and 'shrink' methods. The detone correlation matrix is singular, so it
        cannot be inverted.
    mkt_comp : int, optional
        Number of first components that will be removed using the detone method.

    Returns
    -------
    kurt : nd-array
        The estimation of cokurtosis square matrix.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    """

    if not isinstance(X, pd.DataFrame):
        raise ValueError("X must be a DataFrame")

    assets = X.columns.tolist()
    cols = list(product(assets, assets))
    cols = [str(y) + " - " + str(x) for x, y in cols]

    if method == "hist":
        kurt = cf.cokurtosis_matrix(X)
    if method == "semi":
        kurt = cf.semi_cokurtosis_matrix(X)
    elif method in ["fixed", "spectral", "shrink"]:
        kurt = cf.cokurtosis_matrix(X)
        T, N = X.shape
        q = T / N
        kurt = af.denoiseCov(
            kurt,
            q,
            kind=method,
            bWidth=bWidth,
            detone=detone,
            mkt_comp=mkt_comp,
            alpha=alpha,
        )

    kurt = pd.DataFrame(np.array(kurt, ndmin=2), columns=cols, index=cols)

    return kurt


def forward_regression(X, y, criterion="pvalue", threshold=0.05, verbose=False):
    r"""
    Select the variables that estimate the best model using stepwise
    forward regression. In case none of the variables has a p-value lower
    than threshold, the algorithm will select the variable with lowest p-value.

    Parameters
    ----------
    X : DataFrame of shape (n_samples, n_factors)
        Risk factors returns matrix, where n_samples is the number of samples
        and n_factors is the number of risk factors.
    y : Series of shape (n_samples, 1)
        Asset returns column DataFrame or Series, where n_samples is the number
        of samples.
    criterion : str, optional
        The default is 'pvalue'. Possible values of the criterion used to select
        the best features are:

        - 'pvalue': select the features based on p-values.
        - 'AIC': select the features based on lowest Akaike Information Criterion.
        - 'SIC': select the features based on lowest Schwarz Information Criterion.
        - 'R2': select the features based on highest R Squared.
        - 'R2_A': select the features based on highest Adjusted R Squared.

    threshold : scalar, optional
        Is the maximum p-value for each variable that will be
        accepted in the model. The default is 0.05.
    verbose : bool, optional
        Enable verbose output. The default is False.

    Returns
    -------
    value : list
        A list of the variables that produce the best model.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    """
    if not isinstance(X, pd.DataFrame):
        raise ValueError("X must be a DataFrame")

    if not isinstance(y, pd.DataFrame) and not isinstance(y, pd.Series):
        raise ValueError("y must be a column DataFrame")

    if isinstance(y, pd.DataFrame):
        if y.shape[0] > 1 and y.shape[1] > 1:
            raise ValueError("y must be a column DataFrame")

    included = []
    aic = 1e10
    sic = 1e10
    r2 = -1e10
    r2_a = -1e10
    pvalues = None

    if criterion == "pvalue":
        value = 0
        while value <= threshold:
            excluded = list(set(X.columns) - set(included))
            best_pvalue = 999999
            new_feature = None
            for i in excluded:
                factors = included + [i]
                X1 = X[factors]
                X1 = sm.add_constant(X1)
                results = sm.OLS(y, X1).fit()
                new_pvalues = results.pvalues
                new_pvalues = new_pvalues[new_pvalues.index != "const"]
                cond_1 = new_pvalues.max()
                if best_pvalue > new_pvalues[i] and cond_1 <= threshold:
                    best_pvalue = results.pvalues[i]
                    new_feature = i
                    pvalues = new_pvalues.copy()

            if pvalues is not None:
                value = pvalues[pvalues.index != "const"].max()

            if new_feature is None:
                break
            else:
                included.append(new_feature)

            if verbose:
                print("Add {} with p-value {:.6}".format(new_feature, best_pvalue))

        # This part is how to deal when there isn't an asset with pvalue lower than threshold
        if len(included) == 0:
            excluded = list(set(X.columns) - set(included))
            best_pvalue = 999999
            new_feature = None
            for i in excluded:
                factors = included + [i]
                X1 = X[factors]
                X1 = sm.add_constant(X1)
                results = sm.OLS(y, X1).fit()
                new_pvalues = results.pvalues
                new_pvalues = new_pvalues[new_pvalues.index != "const"]
                if best_pvalue > new_pvalues[i]:
                    best_pvalue = results.pvalues[i]
                    new_feature = i
                    pvalues = new_pvalues.copy()

            value = pvalues[pvalues.index != "const"].max()

            included.append(new_feature)

            if verbose:
                print(
                    "Add {} with p-value {:.6}".format(pvalues.idxmax(), pvalues.max())
                )

    else:
        excluded = X.columns.tolist()
        flag = False
        n = len(excluded)

        for j in range(n):
            value = {}
            n_ini = len(excluded)
            for i in excluded:
                factors = included.copy()
                factors.append(i)
                X1 = X[factors]
                X1 = sm.add_constant(X1)
                results = sm.OLS(y, X1).fit()

                if criterion == "AIC":
                    value[i] = results.aic
                elif criterion == "SIC":
                    value[i] = results.bic
                elif criterion == "R2":
                    value[i] = results.rsquared
                elif criterion == "R2_A":
                    value[i] = results.rsquared_adj

            value = pd.Series(value)

            if criterion in ["AIC", "SIC"]:
                key = value.idxmin()
                value = value.min()
            if criterion in ["R2", "R2_A"]:
                key = value.idxmax()
                value = value.max()

            if criterion == "AIC":
                if value < aic:
                    excluded.remove(key)
                    included.append(key)
                    aic = value
                    flag = True
            elif criterion == "SIC":
                if value < sic:
                    excluded.remove(key)
                    included.append(key)
                    sic = value
                    flag = True
            elif criterion == "R2":
                if value > r2:
                    excluded.remove(key)
                    included.append(key)
                    r2 = value
                    flag = True
            elif criterion == "R2_A":
                if value > r2_a:
                    excluded.remove(key)
                    included.append(key)
                    r2_a = value
                    flag = True

            if n_ini == len(excluded):
                break

            if flag and verbose:
                print("Add {} with {} {:.6}".format(key, criterion, value))

            flag = False

    return included


def backward_regression(X, y, criterion="pvalue", threshold=0.05, verbose=False):
    r"""
    Select the variables that estimate the best model using stepwise
    backward regression. In case none of the variables has a p-value lower
    than threshold, the algorithm will select the variable with lowest p-value.

    Parameters
    ----------
    X : DataFrame of shape (n_samples, n_factors)
        Risk factors returns matrix, where n_samples is the number of samples
        and n_factors is the number of risk factors.
    y : Series of shape (n_samples, 1)
        Asset returns column DataFrame or Series, where n_samples is the number
        of samples.
    criterion : str, optional
        The default is 'pvalue'. Possible values of the criterion used to select
        the best features are:

        - 'pvalue': select the features based on p-values.
        - 'AIC': select the features based on lowest Akaike Information Criterion.
        - 'SIC': select the features based on lowest Schwarz Information Criterion.
        - 'R2': select the features based on highest R Squared.
        - 'R2_A': select the features based on highest Adjusted R Squared.
    threshold : scalar, optional
        Is the maximum p-value for each variable that will be
        accepted in the model. The default is 0.05.
    verbose : bool, optional
        Enable verbose output. The default is False.

    Returns
    -------
    value : list
        A list of the variables that produce the best model.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    """

    if not isinstance(X, pd.DataFrame):
        raise ValueError("X must be a DataFrame")

    if not isinstance(y, pd.DataFrame) and not isinstance(y, pd.Series):
        raise ValueError("y must be a column DataFrame")

    if isinstance(y, pd.DataFrame):
        if y.shape[0] > 1 and y.shape[1] > 1:
            raise ValueError("y must be a column DataFrame")

    X1 = sm.add_constant(X)
    results = sm.OLS(y, X1).fit()
    pvalues = results.pvalues
    aic = results.aic
    sic = results.bic
    r2 = results.rsquared
    r2_a = results.rsquared_adj

    included = pvalues.index.tolist()

    if criterion == "pvalue":
        excluded = ["const"]
        while pvalues[pvalues.index != "const"].max() > threshold:
            factors = pvalues[~pvalues.index.isin(excluded)].index.tolist()
            X1 = X[factors]
            X1 = sm.add_constant(X1)
            results = sm.OLS(y, X1).fit()
            pvalues = results.pvalues
            pvalues = pvalues[pvalues.index != "const"]
            if pvalues.shape[0] == 0:
                break
            excluded = ["const", pvalues.idxmax()]
            if verbose and pvalues.max() > threshold:
                print(
                    "Drop {} with p-value {:.6}".format(pvalues.idxmax(), pvalues.max())
                )

        included = pvalues[pvalues.index != "const"].index.tolist()

        # This part is how to deal when there isn't an asset with pvalue lower than threshold
        if len(included) == 0:
            excluded = list(set(X.columns) - set(included))
            best_pvalue = 999999
            new_feature = None
            for i in excluded:
                factors = included + [i]
                X1 = X[factors]
                X1 = sm.add_constant(X1)
                results = sm.OLS(y, X1).fit()
                new_pvalues = results.pvalues
                new_pvalues = results.pvalues
                new_pvalues = new_pvalues[new_pvalues.index != "const"]
                if best_pvalue > new_pvalues[i]:
                    best_pvalue = results.pvalues[i]
                    new_feature = i
                    pvalues = new_pvalues.copy()

            value = pvalues[pvalues.index != "const"].max()

            included.append(new_feature)

            if verbose:
                print(
                    "Add {} with p-value {:.6}".format(pvalues.idxmax(), pvalues.max())
                )

    else:
        included.remove("const")
        flag = False
        n = len(included)

        for j in range(n):
            value = {}
            n_ini = len(included)
            for i in included:
                factors = included.copy()
                factors.remove(i)
                X1 = X[factors]
                X1 = sm.add_constant(X1)
                results = sm.OLS(y, X1).fit()

                if criterion == "AIC":
                    value[i] = results.aic
                elif criterion == "SIC":
                    value[i] = results.bic
                elif criterion == "R2":
                    value[i] = results.rsquared
                elif criterion == "R2_A":
                    value[i] = results.rsquared_adj

            value = pd.Series(value)

            if criterion in ["AIC", "SIC"]:
                key = value.idxmin()
                value = value.min()
            if criterion in ["R2", "R2_A"]:
                key = value.idxmax()
                value = value.max()

            if criterion == "AIC":
                if value < aic:
                    included.remove(key)
                    aic = value
                    flag = True
            elif criterion == "SIC":
                if value < sic:
                    included.remove(key)
                    sic = value
                    flag = True
            elif criterion == "R2":
                if value > r2:
                    included.remove(key)
                    r2 = value
                    flag = True
            elif criterion == "R2_A":
                if value > r2_a:
                    included.remove(key)
                    r2_a = value
                    flag = True

            if n_ini == len(included):
                break

            if flag and verbose:
                print("Drop {} with {} {:.6}".format(key, criterion, value))

            flag = False

    return included


def PCR(X, y, n_components=0.95):
    r"""
    Estimate the coefficients using Principal Components Regression (PCR).

    Parameters
    ----------
    X : DataFrame of shape (n_samples, n_factors)
        Risk factors returns matrix, where n_samples is the number of samples
        and n_factors is the number of risk factors.
    y : DataFrame or Series of shape (n_samples, 1)
        Asset returns column DataFrame or Series, where n_samples is the number
        of samples.
    n_components : int, float, None or str, optional
        if 1 < n_components (int), it represents the number of components that
        will be keep. if 0 < n_components < 1 (float), it represents the
        percentage of variance that the is explained by the components kept.
        See `PCA <https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.PCA.html>`_
        for more details. The default is 0.95.

    Returns
    -------
    value : nd-array
        An array with the coefficients of the model calculated using PCR.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    """

    if not isinstance(X, pd.DataFrame):
        raise ValueError("X must be a DataFrame")

    if not isinstance(y, pd.DataFrame) and not isinstance(y, pd.Series):
        raise ValueError("y must be a column DataFrame")

    if isinstance(y, pd.DataFrame):
        if y.shape[0] > 1 and y.shape[1] > 1:
            raise ValueError("y must be a column DataFrame")

    scaler = StandardScaler()
    scaler.fit(X)
    X_std = scaler.transform(X)

    if n_components > 0 and n_components < 1:
        pca = PCA(n_components=n_components)
    elif n_components >= 1:
        pca = PCA(n_components=int(n_components))

    pca.fit(X_std)
    Z_p = pca.transform(X_std)
    V_p = pca.components_.T

    results = sm.OLS(y, sm.add_constant(Z_p)).fit()
    beta_pc = results.params[1:]
    beta_pc = np.array(beta_pc, ndmin=2)

    std = np.array(np.std(X, axis=0, ddof=1), ndmin=2)
    mean = np.array(np.mean(X, axis=0), ndmin=2)
    beta = V_p @ beta_pc.T / std.T

    beta_0 = np.array(y.mean(), ndmin=2) - np.sum(beta * mean.T)

    beta = np.insert(beta, 0, beta_0)
    beta = np.array(beta, ndmin=2)

    return beta


def loadings_matrix(
    X,
    Y,
    feature_selection="stepwise",
    stepwise="Forward",
    criterion="pvalue",
    threshold=0.05,
    n_components=0.95,
    verbose=False,
):
    r"""
    Estimate the loadings matrix using stepwise regression.

    Parameters
    ----------
    X : DataFrame of shape (n_samples, n_factors)
        Risk factors returns matrix, where n_samples is the number of samples
        and n_factors is the number of risk factors.
    Y : DataFrame of shape (n_samples, n_assets)
        Assets returns DataFrame, where n_samples is the number of
        observations and n_assets is the number of assets.
    feature_selection: str, 'stepwise' or 'PCR', optional
        Indicate the method used to estimate the loadings matrix.
        The default is 'stepwise'.  Possible values are:

        - 'stepwise': use stepwise regression to select the best factors and estimate coefficients.
        - 'PCR': use principal components regression to estimate coefficients.
    stepwise: str 'Forward' or 'Backward', optional
        Indicate the method used for stepwise regression.
        The default is 'Forward'.
    criterion : str, optional
        The default is 'pvalue'. Possible values of the criterion used to select
        the best features are:

        - 'pvalue': select the features based on p-values.
        - 'AIC': select the features based on lowest Akaike Information Criterion.
        - 'SIC': select the features based on lowest Schwarz Information Criterion.
        - 'R2': select the features based on highest R Squared.
        - 'R2_A': select the features based on highest Adjusted R Squared.
    threshold : scalar, optional
        Is the maximum p-value for each variable that will be
        accepted in the model. The default is 0.05.
    n_components : int, float, None or str, optional
        if 1 < n_components (int), it represents the number of components that
        will be keep. if 0 < n_components < 1 (float), it represents the
        percentage of variance that the is explained by the components kept.
        See `PCA <https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.PCA.html>`_
        for more details. The default is 0.95.
    verbose : bool, optional
        Enable verbose output. The default is False.

    Returns
    -------
    loadings : DataFrame
        Loadings matrix.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    """
    if not isinstance(X, pd.DataFrame):
        raise ValueError("X must be a DataFrame")

    if not isinstance(Y, pd.DataFrame):
        raise ValueError("Y must be a DataFrame")

    rows = Y.columns.tolist()
    cols = X.columns.tolist()
    cols.insert(0, "const")
    loadings = np.zeros((len(rows), len(cols)))
    loadings = pd.DataFrame(loadings, index=rows, columns=cols)

    for i in rows:
        if feature_selection == "stepwise":
            if stepwise == "Forward":
                included = forward_regression(
                    X, Y[i], criterion=criterion, threshold=threshold, verbose=verbose
                )
            elif stepwise == "Backward":
                included = backward_regression(
                    X, Y[i], criterion=criterion, threshold=threshold, verbose=verbose
                )
            else:
                raise ValueError("Choose and adequate stepwise method")
            results = sm.OLS(Y[i], sm.add_constant(X[included])).fit()
            params = results.params
            loadings.loc[i, params.index.tolist()] = params.T
        elif feature_selection == "PCR":
            beta = PCR(X, Y[i], n_components=n_components)
            beta = pd.Series(np.ravel(beta), index=cols)
            loadings.loc[i, cols] = beta.T

    return loadings


def risk_factors(
    X,
    Y,
    B=None,
    const=True,
    method_mu="hist",
    method_cov="hist",
    feature_selection="stepwise",
    stepwise="Forward",
    criterion="pvalue",
    threshold=0.05,
    n_components=0.95,
    dict_mu={},
    dict_cov={},
):
    r"""
    Estimate the expected returns vector and covariance matrix based on risk
    factors models :cite:`b-Ross` :cite:`b-Fan`.

    .. math::
        \begin{aligned}
        R & = \alpha + B F + \epsilon \\
        \mu_{f} & = \alpha +BE(F) \\
        \Sigma_{f} & = B \Sigma_{F} B^{T} + \Sigma_{\epsilon} \\
        \end{aligned}


    where:

    :math:`R` is the series returns.

    :math:`\alpha` is the intercept.

    :math:`B` is the loadings matrix.

    :math:`F` is the expected returns vector of the risk factors.

    :math:`\Sigma_{F}` is the covariance matrix of the risk factors.

    :math:`\Sigma_{\epsilon}` is the covariance matrix of error terms.

    :math:`\mu_{f}` is the expected returns vector obtained with the
    risk factor model.

    :math:`\Sigma_{f}` is the covariance matrix obtained with the risk
    factor model.

    Parameters
    ----------
    X : DataFrame of shape (n_samples, n_factors)
        Risk factors returns matrix, where n_samples is the number of samples
        and n_factors is the number of risk factors.
    Y : DataFrame of shape (n_samples, n_assets)
        Assets returns DataFrame, where n_samples is the number of
        observations and n_assets is the number of assets.
    B : DataFrame of shape (n_assets, n_factors), optional
        Loadings matrix, where n_assets is the number assets and n_factors is
        the number of risk factors. If is not specified, is estimated using
        stepwise regression. The default is None.
    const : bool, optional
        Indicate if the loadings matrix has a constant.
        The default is False.
    method_mu : str, optional
        The method used to estimate the expected returns of factors.
        The default value is 'hist'. Possible values are:

        - 'hist': use historical estimates.
        - 'ewma1'': use ewma with adjust=True, see `EWM <https://pandas.pydata.org/pandas-docs/stable/user_guide/computation.html#exponentially-weighted-windows>`__ for more details.
        - 'ewma2': use ewma with adjust=False, see `EWM <https://pandas.pydata.org/pandas-docs/stable/user_guide/computation.html#exponentially-weighted-windows>`__ for more details.
        - 'JS': James-Stein estimator. For more information see :cite:`b-Meucci2005` and :cite:`b-Feng2016`.
        - 'BS': Bayes-Stein estimator. For more information see :cite:`b-Jorion1986`.
        - 'BOP': BOP estimator. For more information see :cite:`b-Bodnar2019`.
    method_cov : str, optional
        The method used to estimate the covariance matrix of factors.
        The default is 'hist'. Possible values are:

        - 'hist': use historical estimates.
        - 'ewma1'': use ewma with adjust=True, see `EWM <https://pandas.pydata.org/pandas-docs/stable/user_guide/computation.html#exponentially-weighted-windows>`__ for more details.
        - 'ewma2': use ewma with adjust=False, see `EWM <https://pandas.pydata.org/pandas-docs/stable/user_guide/computation.html#exponentially-weighted-windows>`__ for more details.
        - 'ledoit': use the Ledoit and Wolf Shrinkage method.
        - 'oas': use the Oracle Approximation Shrinkage method.
        - 'shrunk': use the basic Shrunk Covariance method.
        - 'gl': use the basic Graphical Lasso Covariance method.
        - 'jlogo': use the j-LoGo Covariance method. For more information see: :cite:`b-jLogo`.
        - 'fixed': denoise using fixed method. For more information see chapter 2 of :cite:`b-MLforAM`.
        - 'spectral': denoise using spectral method. For more information see chapter 2 of :cite:`b-MLforAM`.
        - 'shrink': denoise using shrink method. For more information see chapter 2 of :cite:`b-MLforAM`.
        - 'gerber1': use the Gerber statistic 1. For more information see: :cite:`b-Gerber2021`.
        - 'gerber2': use the Gerber statistic 2. For more information see: :cite:`b-Gerber2021`.
    feature_selection: str, 'stepwise' or 'PCR', optional
        Indicate the method used to estimate the loadings matrix.
        The default is 'stepwise'.  Possible values are:

        - 'stepwise': use stepwise regression to select the best factors and estimate coefficients.
        - 'PCR': use principal components regression to estimate coefficients.
    stepwise: str, 'Forward' or 'Backward'
        Indicate the method used for stepwise regression.
        The default is 'Forward'.
    criterion : str, optional
        The default is 'pvalue'. Possible values of the criterion used to select
        the best features are:

        - 'pvalue': select the features based on p-values.
        - 'AIC': select the features based on lowest Akaike Information Criterion.
        - 'SIC': select the features based on lowest Schwarz Information Criterion.
        - 'R2': select the features based on highest R Squared.
        - 'R2_A': select the features based on highest Adjusted R Squared.
    threshold : scalar, optional
        Is the maximum p-value for each variable that will be
        accepted in the model. The default is 0.05.
    n_components : int, float, None or str, optional
        if 1 < n_components (int), it represents the number of components that
        will be keep. if 0 < n_components < 1 (float), it represents the
        percentage of variance that the is explained by the components kept.
        See `PCA <https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.PCA.html>`_
        for more details. The default is 0.95.
    dict_mu : dict
        Other variables related to the expected returns.
    dict_cov : dict
        Other variables related to the covariance estimation.

    Returns
    -------
    mu : DataFrame
        The mean vector of risk factors model.
    cov : DataFrame
        The covariance matrix of risk factors model.
    returns : DataFrame
        The returns based on a risk factor model.
    B : DataFrame
        Loadings matrix.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    """
    if not isinstance(X, pd.DataFrame) and not isinstance(Y, pd.DataFrame):
        raise ValueError("X and Y must be DataFrames")

    if B is None:
        B = loadings_matrix(
            X,
            Y,
            feature_selection=feature_selection,
            stepwise=stepwise,
            criterion=criterion,
            threshold=threshold,
            n_components=n_components,
            verbose=False,
        )
    elif not isinstance(B, pd.DataFrame):
        raise ValueError("B must be a DataFrame")

    assets = Y.columns.tolist()
    dates = X.index.tolist()

    X1 = X.copy()
    if const == True or ("const" in B.columns.tolist()):
        mu_f = np.hstack(
            [
                np.ones((1, 1)),
                np.array(mean_vector(X1, method=method_mu, **dict_mu), ndmin=2),
            ]
        )
        X1 = sm.add_constant(X)
    else:
        mu_f = np.array(mean_vector(X1, method=method_mu, **dict_mu), ndmin=2)
    S_f = np.array(covar_matrix(X1, method=method_cov, **dict_cov), ndmin=2)
    B_ = np.array(B, ndmin=2)

    returns = np.array(X1, ndmin=2) @ B_.T
    mu = B_ @ mu_f.T

    e = np.array(Y, ndmin=2) - returns
    S_e = np.diag(np.var(np.array(e), ddof=1, axis=0))
    S = B_ @ S_f @ B_.T + S_e

    mu = pd.DataFrame(mu.T, columns=assets)
    cov = pd.DataFrame(S, index=assets, columns=assets)
    returns = pd.DataFrame(returns, index=dates, columns=assets)

    return mu, cov, returns, B


def black_litterman(
    X,
    w,
    P,
    Q,
    delta=1,
    rf=0,
    eq=True,
    method_mu="hist",
    method_cov="hist",
    dict_mu={},
    dict_cov={},
):
    r"""
    Estimate the expected returns vector and covariance matrix based
    on the Black Litterman model :cite:`b-BlackLitterman` :cite:`b-Black1`.

    .. math::
        \begin{aligned}
        \Pi & = \delta \Sigma w \\
        \Pi_{BL} & = \left [ (\tau\Sigma)^{-1}+ P^{T} \Omega^{-1}P \right]^{-1}
        \left[(\tau\Sigma)^{-1} \Pi + P^{T} \Omega^{-1} Q \right] \\
        M & = \left((\tau\Sigma)^{-1} + P^{T}\Omega^{-1} P \right)^{-1} \\
        \mu_{BL} & = \Pi_{BL} + r_{f} \\
        \Sigma_{BL} & = \Sigma + M \\
        \end{aligned}


    where:

    :math:`r_{f}` is the risk free rate.

    :math:`\delta` is the risk aversion factor.

    :math:`\Pi` is the equilibrium excess returns.

    :math:`\Sigma` is the covariance matrix.

    :math:`P` is the views matrix.

    :math:`Q` is the views returns matrix.

    :math:`\Omega` is the covariance matrix of the error views.

    :math:`\mu_{BL}` is the mean vector obtained with the black
    litterman model.

    :math:`\Sigma_{BL}` is the covariance matrix obtained with the black
    litterman model.

    Parameters
    ----------
    X : DataFrame of shape (n_samples, n_assets)
        Assets returns DataFrame, where n_samples is the number of
        observations and n_assets is the number of assets.
    w : DataFrame or Series of shape (n_assets, 1)
        Portfolio weights, where n_assets is the number of assets.
    P : DataFrame of shape (n_views, n_assets)
        Analyst's views matrix, can be relative or absolute.
    Q : DataFrame of shape (n_views, 1)
        Expected returns of analyst's views.
    delta : float, optional
        Risk aversion factor. The default value is 1.
    rf : scalar, optional
        Risk free rate. The default is 0.
    eq : bool, optional
        Indicate if use equilibrium or historical excess returns.
        The default is True.
    method_mu : str, optional
        The method used to estimate the expected returns.
        The default value is 'hist'.

        - 'hist': use historical estimates.
        - 'ewma1': use ewma with adjust=True. For more information see `EWM <https://pandas.pydata.org/pandas-docs/stable/user_guide/window.html#exponentially-weighted-window>`__.
        - 'ewma2': use ewma with adjust=False. For more information see `EWM <https://pandas.pydata.org/pandas-docs/stable/user_guide/window.html#exponentially-weighted-window>`__.
        - 'JS': James-Stein estimator. For more information see :cite:`b-Meucci2005` and :cite:`b-Feng2016`.
        - 'BS': Bayes-Stein estimator. For more information see :cite:`b-Jorion1986`.
        - 'BOP': BOP estimator. For more information see :cite:`b-Bodnar2019`.
    method_cov : str, optional
        The method used to estimate the covariance matrix.
        The default is 'hist'. Possible values are:

        - 'hist': use historical estimates.
        - 'ewma1': use ewma with adjust=True. For more information see `EWM <https://pandas.pydata.org/pandas-docs/stable/user_guide/window.html#exponentially-weighted-window>`__.
        - 'ewma2': use ewma with adjust=False. For more information see `EWM <https://pandas.pydata.org/pandas-docs/stable/user_guide/window.html#exponentially-weighted-window>`__.
        - 'ledoit': use the Ledoit and Wolf Shrinkage method.
        - 'oas': use the Oracle Approximation Shrinkage method.
        - 'shrunk': use the basic Shrunk Covariance method.
        - 'gl': use the basic Graphical Lasso Covariance method.
        - 'jlogo': use the j-LoGo Covariance method. For more information see: :cite:`b-jLogo`.
        - 'fixed': denoise using fixed method. For more information see chapter 2 of :cite:`b-MLforAM`.
        - 'spectral': denoise using spectral method. For more information see chapter 2 of :cite:`b-MLforAM`.
        - 'shrink': denoise using shrink method. For more information see chapter 2 of :cite:`b-MLforAM`.
        - 'gerber1': use the Gerber statistic 1. For more information see: :cite:`b-Gerber2021`.
        - 'gerber2': use the Gerber statistic 2. For more information see: :cite:`b-Gerber2021`.
    dict_mu : dict
        Other variables related to the mean vector estimation method.
    dict_cov : dict
        Other variables related to the covariance estimation method.

    Returns
    -------
    mu : DataFrame
        The mean vector of Black Litterman model.
    cov : DataFrame
        The covariance matrix of Black Litterman model.
    w : DataFrame
        The equilibrium weights of Black Litterman model, without constraints.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    """
    if not isinstance(X, pd.DataFrame) and not isinstance(w, pd.DataFrame):
        raise ValueError("X and w must be DataFrames")

    if w.shape[0] > 1 and w.shape[1] > 1:
        raise ValueError("w must be a column DataFrame")

    assets = X.columns.tolist()

    w = np.array(w, ndmin=2)
    if w.shape[0] == 1:
        w = w.T

    mu = np.array(mean_vector(X, method=method_mu, **dict_mu), ndmin=2)
    S = np.array(covar_matrix(X, method=method_cov, **dict_cov), ndmin=2)
    P = np.array(P, ndmin=2)
    Q = np.array(Q, ndmin=2)
    tau = 1 / X.shape[0]
    Omega = np.array(np.diag(np.diag(P @ (tau * S) @ P.T)), ndmin=2)

    if eq == True:
        PI = delta * (S @ w)
    elif eq == False:
        PI = mu.T - rf

    PI_ = inv(inv(tau * S) + P.T @ inv(Omega) @ P) @ (
        inv(tau * S) @ PI + P.T @ inv(Omega) @ Q
    )
    M = inv(inv(tau * S) + P.T @ inv(Omega) @ P)
    # PI_1 = PI + (tau * S* P.T) * inv(P * tau * S * P.T + Omega) * (Q - P * PI)
    # M = tau * S - (tau * S * P.T) * inv(P * tau * S * P.T + Omega) * P * tau * S

    mu = PI_ + rf
    mu = mu.T
    cov = S + M
    w = inv(delta * cov) @ PI_

    mu = pd.DataFrame(mu, columns=assets)
    cov = pd.DataFrame(cov, index=assets, columns=assets)
    w = pd.DataFrame(w, index=assets)

    return mu, cov, w


def augmented_black_litterman(
    X,
    w,
    F,
    B,
    P=None,
    Q=None,
    P_f=None,
    Q_f=None,
    delta=1,
    rf=0,
    eq=True,
    const=True,
    method_mu="hist",
    method_cov="hist",
    dict_mu={},
    dict_cov={},
):
    r"""
    Estimate the expected returns vector and covariance matrix based
    on the Augmented Black Litterman model :cite:`b-WCheung`.

    .. math::
        \begin{aligned}
        \Pi^{a} & = \delta \left [ \begin{array}{c} \Sigma \\ \Sigma_{F} B^{T} \\ \end{array} \right ] w \\
        P^{a} & = \left [ \begin{array}{cc} P & 0 \\ 0 & P_{F} \\ \end{array} \right ] \\
        Q^{a} & = \left [ \begin{array}{c} Q \\ Q_{F} \\ \end{array} \right ] \\
        \Sigma^{a} & = \left [ \begin{array}{cc} \Sigma & B \Sigma_{F}\\ \Sigma_{F} B^{T} & \Sigma_{F} \\ \end{array} \right ] \\
        \Omega^{a} & = \left [ \begin{array}{cc} \Omega & 0 \\ 0 & \Omega_{F} \\ \end{array} \right ] \\
        \Pi^{a}_{BL} & = \left [ (\tau \Sigma^{a})^{-1} + (P^{a})^{T} (\Omega^{a})^{-1} P^{a} \right ]^{-1}
        \left [ (\tau\Sigma^{a})^{-1} \Pi^{a} + (P^{a})^{T} (\Omega^{a})^{-1} Q^{a} \right ] \\
        M^{a} & = \left ( (\tau\Sigma^{a})^{-1} + (P^{a})^{T} (\Omega^{a})^{-1} P^{a} \right )^{-1} \\
        \mu^{a}_{BL} & = \Pi^{a}_{BL} + r_{f} \\
        \Sigma^{a}_{BL} & = \Sigma^{a} + M^{a} \\
        \end{aligned}


    where:

    :math:`r_{f}` is the risk free rate.

    :math:`\delta` is the risk aversion factor.

    :math:`B` is the loadings matrix.

    :math:`\Sigma` is the covariance matrix of assets.

    :math:`\Sigma_{F}` is the covariance matrix of factors.

    :math:`\Sigma^{a}` is the augmented covariance matrix.

    :math:`P` is the assets views matrix.

    :math:`Q` is the assets views returns matrix.

    :math:`P_{F}` is the factors views matrix.

    :math:`Q_{F}` is the factors views returns matrix.

    :math:`P^{a}` is the augmented views matrix.

    :math:`Q^{a}` is the augmented views returns matrix.

    :math:`\Pi^{a}` is the augmented equilibrium excess returns.

    :math:`\Omega` is the covariance matrix of errors of assets views.

    :math:`\Omega_{F}` is the covariance matrix of errors of factors views.

    :math:`\Omega^{a}` is the covariance matrix of errors of augmented views.

    :math:`\mu^{a}_{BL}` is the mean vector obtained with the Augmented Black
    Litterman model.

    :math:`\Sigma^{a}_{BL}` is the covariance matrix obtained with the Augmented
    Black Litterman model.

    Parameters
    ----------
    X : DataFrame of shape (n_samples, n_assets)
        Assets returns DataFrame, where n_samples is the number of
        observations and n_assets is the number of assets.
    w : DataFrame or Series of shape (n_assets, 1)
        Portfolio weights, where n_assets is the number of assets.
    F : DataFrame of shape (n_samples, n_factors)
        Risk factors returns DataFrame, where n_samples is the number of samples
        and n_factors is the number of risk factors.
    B : DataFrame of shape (n_assets, n_factors), optional
        Loadings matrix, where n_assets is the number assets and n_factors is
        the number of risk factors.
    P : DataFrame of shape (n_views, n_assets)
        Analyst's views matrix, can be relative or absolute.
    Q : DataFrame of shape (n_views, 1)
        Expected returns of analyst's views.
    P_f : DataFrame of shape (n_views, n_factors)
        Analyst's factors views matrix, can be relative or absolute.
    Q_f : DataFrame of shape (n_views, 1)
        Expected returns of analyst's factors views.
    delta : float, optional
        Risk aversion factor. The default value is 1.
    rf : scalar, optional
        Risk free rate. The default is 0.
    eq : bool, optional
        Indicate if use equilibrium or historical excess returns.
        The default is True.
    const : bool, optional
        Indicate if the loadings matrix has a constant.
        The default is True.
    method_mu : str, optional
        The method used to estimate the expected returns.
        The default value is 'hist'.

        - 'hist': use historical estimates.
        - 'ewma1': use ewma with adjust=True. For more information see `EWM <https://pandas.pydata.org/pandas-docs/stable/user_guide/window.html#exponentially-weighted-window>`__.
        - 'ewma2': use ewma with adjust=False. For more information see `EWM <https://pandas.pydata.org/pandas-docs/stable/user_guide/window.html#exponentially-weighted-window>`__.
        - 'JS': James-Stein estimator. For more information see :cite:`b-Meucci2005` and :cite:`b-Feng2016`.
        - 'BS': Bayes-Stein estimator. For more information see :cite:`b-Jorion1986`.
        - 'BOP': BOP estimator. For more information see :cite:`b-Bodnar2019`.
    method_cov : str, optional
        The method used to estimate the covariance matrix.
        The default is 'hist'. Possible values are:

        - 'hist': use historical estimates.
        - 'ewma1': use ewma with adjust=True. For more information see `EWM <https://pandas.pydata.org/pandas-docs/stable/user_guide/window.html#exponentially-weighted-window>`__.
        - 'ewma2': use ewma with adjust=False. For more information see `EWM <https://pandas.pydata.org/pandas-docs/stable/user_guide/window.html#exponentially-weighted-window>`__.
        - 'ledoit': use the Ledoit and Wolf Shrinkage method.
        - 'oas': use the Oracle Approximation Shrinkage method.
        - 'shrunk': use the basic Shrunk Covariance method.
        - 'gl': use the basic Graphical Lasso Covariance method.
        - 'jlogo': use the j-LoGo Covariance method. For more information see: :cite:`b-jLogo`.
        - 'fixed': denoise using fixed method. For more information see chapter 2 of :cite:`b-MLforAM`.
        - 'spectral': denoise using spectral method. For more information see chapter 2 of :cite:`b-MLforAM`.
        - 'shrink': denoise using shrink method. For more information see chapter 2 of :cite:`b-MLforAM`.
        - 'gerber1': use the Gerber statistic 1. For more information see: :cite:`b-Gerber2021`.
        - 'gerber2': use the Gerber statistic 2. For more information see: :cite:`b-Gerber2021`.
    dict_mu : dict
        Other variables related to the mean vector estimation method.
    dict_cov : dict
        Other variables related to the covariance estimation method.

    Returns
    -------
    mu : DataFrame
        The mean vector of Augmented Black Litterman model.
    cov : DataFrame
        The covariance matrix of Augmented Black Litterman model.
    w : DataFrame
        The equilibrium weights of Augmented Black Litterman model, without constraints.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    """
    if not isinstance(X, pd.DataFrame) and not isinstance(w, pd.DataFrame):
        raise ValueError("X and w must be DataFrames")

    if not isinstance(F, pd.DataFrame) and not isinstance(B, pd.DataFrame):
        raise ValueError("F and B must be DataFrames")

    if w.shape[0] > 1 and w.shape[1] > 1:
        raise ValueError("w must be a column DataFrame")

    assets = X.columns.tolist()
    N = len(assets)

    w = np.array(w, ndmin=2)
    if w.shape[0] == 1:
        w = w.T

    if B is not None:
        B_ = np.array(B, ndmin=2)
        if const == True:
            alpha = B_[:, :1]
            B_ = B_[:, 1:]

    mu = np.array(mean_vector(X, method=method_mu, **dict_mu), ndmin=2)
    S = np.array(covar_matrix(X, method=method_cov, **dict_cov), ndmin=2)

    tau = 1 / X.shape[0]

    mu_f = np.array(mean_vector(F, method=method_mu, **dict_mu), ndmin=2)
    S_f = np.array(covar_matrix(F, method=method_cov, **dict_cov), ndmin=2)

    if P is not None and Q is not None and P_f is None and Q_f is None:
        S_a = S
        P_a = P
        Q_a = Q
        Omega = np.array(np.diag(np.diag(P @ (tau * S) @ P.T)), ndmin=2)
        Omega_a = Omega

        if eq == True:
            PI_a_ = delta * S_a @ w
        elif eq == False:
            PI_a_ = mu.T - rf
    elif P is None and Q is None and P_f is not None and Q_f is not None:
        S_a = S_f
        P_a = P_f
        Q_a = Q_f
        Omega_f = np.array(np.diag(np.diag(P_f @ (tau * S_f) @ P_f.T)), ndmin=2)
        Omega_a = Omega_f

        if eq == True:
            PI_a_ = delta * (S_f @ B.T) @ w
        elif eq == False:
            PI_a_ = mu_f.T - rf

    elif P is not None and Q is not None and P_f is not None and Q_f is not None:
        S_a = np.hstack((np.vstack((S, S_f @ B_.T)), np.vstack((B_ @ S_f, S_f))))

        P = np.array(P, ndmin=2)
        Q = np.array(Q, ndmin=2)
        P_f = np.array(P_f, ndmin=2)
        Q_f = np.array(Q_f, ndmin=2)
        zeros_1 = np.zeros((P_f.shape[0], P.shape[1]))
        zeros_2 = np.zeros((P.shape[0], P_f.shape[1]))
        P_a = np.hstack((np.vstack((P, zeros_1)), np.vstack((zeros_2, P_f))))
        Q_a = np.vstack((Q, Q_f))

        Omega = np.array(np.diag(np.diag(P @ (tau * S) @ P.T)), ndmin=2)
        Omega_f = np.array(np.diag(np.diag(P_f @ (tau * S_f) @ P_f.T)), ndmin=2)
        zeros = np.zeros((Omega.shape[0], Omega_f.shape[0]))
        Omega_a = np.hstack((np.vstack((Omega, zeros.T)), np.vstack((zeros, Omega_f))))

        if eq == True:
            PI_a_ = delta * (np.vstack((S, S_f @ B_.T)) @ w)
        elif eq == False:
            PI_a_ = np.vstack((mu.T, mu_f.T)) - rf

    PI_a = inv(inv(tau * S_a) + P_a.T @ inv(Omega_a) @ P_a) @ (
        inv(tau * S_a) @ PI_a_ + P_a.T @ inv(Omega_a) @ Q_a
    )
    M_a = inv(inv(tau * S_a) + P_a.T @ inv(Omega_a) @ P_a)
    # PI_a = PI_a_ + (tau * S_a @ P_a.T) * inv(P_a @ tau * S_a @ P_a.T + Omega) * (Q_a - P_a @ PI_a_)
    # M = tau * S_a - (tau * S_a @ P_a.T) * inv(P_a @ tau * S_a @ P_a.T + Omega_a) @ P_a @ tau * S_a

    mu_a = PI_a + rf
    mu_a = mu_a.T
    cov_a = S_a + M_a
    w_a = inv(delta * cov_a) @ PI_a

    if P is None and Q is None and P_f is not None and Q_f is not None:
        mu_a = mu_a @ B_.T
        cov_a = B_ @ cov_a @ B_.T
        w_a = inv(delta * cov_a) @ B_ @ PI_a

    if const == True:
        mu_a = mu_a[:, :N] + alpha.T

    mu_a = pd.DataFrame(mu_a[:, :N], columns=assets)
    cov_a = pd.DataFrame(cov_a[:N, :N], index=assets, columns=assets)
    w_a = pd.DataFrame(w_a[:N, 0], index=assets)

    return mu_a, cov_a, w_a


def black_litterman_bayesian(
    X,
    F,
    B,
    P_f,
    Q_f,
    delta=1,
    rf=0,
    eq=True,
    const=True,
    method_mu="hist",
    method_cov="hist",
    dict_mu={},
    dict_cov={},
):
    r"""
    Estimate the expected returns vector and covariance matrix based
    on the black litterman model :cite:`b-BLB`.

    .. math::
        \begin{aligned}
        \Sigma_{F} & = B \Sigma_{F} B^{T} + D \\
        \overline{\Pi}_{F} & = \left ( \Sigma_{F}^{-1} + P_{F}^{T}\Omega_{F}^{-1}P_{F} \right )^{-1} \left ( \Sigma_{F}^{-1}\Pi_{F} + P_{F}^{T}\Omega_{F}^{-1}Q_{F} \right) \\
        \overline{\Sigma}_{F} & = \left ( \Sigma_{F}^{-1} + P_{F}^{T}\Omega_{F}^{-1}P_{F} \right )^{-1} \\
        \Sigma_{BLB} & = \left( \Sigma^{-1} - \Sigma^{-1} B \left( \overline{\Sigma}_{F}^{-1} + B^{T}\Sigma^{-1}B \right)^{-1} B^{T}\Sigma^{-1} \right )^{-1} \\
        \mu_{BLB} & = \Sigma_{BLB} \left ( \Sigma^{-1} B \left( \overline{\Sigma}_{F}^{-1} +B^{T}\Sigma^{-1}B \right)^{-1} \overline{\Sigma}_{F}^{-1} \overline{\Pi}_{F} \right ) + r_{f} \\
        \end{aligned}


    where:

    :math:`r_{f}` is the risk free rate.

    :math:`B` is the loadings matrix.

    :math:`D` is a diagonal matrix of variance of errors of a factor model.

    :math:`\Sigma` is the covariance matrix obtained with a factor model.

    :math:`\Pi_{F}` is the equilibrium excess returns of factors.

    :math:`\overline{\Pi}_{F}` is the posterior excess returns of factors.

    :math:`\Sigma_{F}` is the covariance matrix of factors.

    :math:`\overline{\Sigma}_{F}` is the posterior covariance matrix of factors.

    :math:`P_{F}` is the factors views matrix.

    :math:`Q_{F}` is the factors views returns matrix.

    :math:`\Omega_{F}` is the covariance matrix of errors of factors views.

    :math:`\mu_{BLB}` is the mean vector obtained with the Black
    Litterman Bayesian model or posterior predictive mean.

    :math:`\Sigma_{BLB}` is the covariance matrix obtained with the Black
    Litterman Bayesian model or posterior predictive covariance.

    Parameters
    ----------
    X : DataFrame of shape (n_samples, n_assets)
        Assets returns DataFrame, where n_samples is the number of
        observations and n_assets is the number of assets.
    F : DataFrame of shape (n_samples, n_factors)
        Risk factors returns DataFrame, where n_samples is the number of samples
        and n_factors is the number of risk factors.
    B : DataFrame of shape (n_assets, n_factors), optional
        Loadings matrix, where n_assets is the number assets and n_factors is
        the number of risk factors. The default is None.
    P_f : DataFrame of shape (n_views, n_factors)
        Analyst's factors views matrix, can be relative or absolute.
    Q_f : DataFrame of shape (n_views, 1)
        Expected returns of analyst's factors views.
    delta : float, optional
        Risk aversion factor. The default value is 1.
    rf : scalar, optional
        Risk free rate. The default is 0.
    eq : bool, optional
        Indicate if use equilibrium or historical excess returns.
        The default is True.
    const : bool, optional
        Indicate if the loadings matrix has a constant.
        The default is True.
    method_mu : str, optional
        The method used to estimate the expected returns.
        The default value is 'hist'.

        - 'hist': use historical estimates.
        - 'ewma1': use ewma with adjust=True. For more information see `EWM <https://pandas.pydata.org/pandas-docs/stable/user_guide/window.html#exponentially-weighted-window>`__.
        - 'ewma2': use ewma with adjust=False, For more information see `EWM <https://pandas.pydata.org/pandas-docs/stable/user_guide/window.html#exponentially-weighted-window>`__.
        - 'JS': James-Stein estimator. For more information see :cite:`b-Meucci2005` and :cite:`b-Feng2016`.
        - 'BS': Bayes-Stein estimator. For more information see :cite:`b-Jorion1986`.
        - 'BOP': BOP estimator. For more information see :cite:`b-Bodnar2019`.
    method_cov : str, optional
        The method used to estimate the covariance matrix:
        The default is 'hist'. Possible values are:

        - 'hist': use historical estimates.
        - 'ewma1': use ewma with adjust=True. For more information see `EWM <https://pandas.pydata.org/pandas-docs/stable/user_guide/window.html#exponentially-weighted-window>`__.
        - 'ewma2': use ewma with adjust=False. For more information see `EWM <https://pandas.pydata.org/pandas-docs/stable/user_guide/window.html#exponentially-weighted-window>`__.
        - 'ledoit': use the Ledoit and Wolf Shrinkage method.
        - 'oas': use the Oracle Approximation Shrinkage method.
        - 'shrunk': use the basic Shrunk Covariance method.
        - 'gl': use the basic Graphical Lasso Covariance method.
        - 'jlogo': use the j-LoGo Covariance method. For more information see: :cite:`b-jLogo`.
        - 'fixed': denoise using fixed method. For more information see chapter 2 of :cite:`b-MLforAM`.
        - 'spectral': denoise using spectral method. For more information see chapter 2 of :cite:`b-MLforAM`.
        - 'shrink': denoise using shrink method. For more information see chapter 2 of :cite:`b-MLforAM`.
        - 'gerber1': use the Gerber statistic 1. For more information see: :cite:`b-Gerber2021`.
        - 'gerber2': use the Gerber statistic 2. For more information see: :cite:`b-Gerber2021`.
    dict_mu : dict
        Other variables related to the mean vector estimation method.
    dict_cov : dict
        Other variables related to the covariance estimation method.

    Returns
    -------
    mu : DataFrame
        The mean vector of Black Litterman model.
    cov : DataFrame
        The covariance matrix of Black Litterman model.
    w : DataFrame
        The equilibrium weights of Black Litterman model, without constraints.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    """
    if not isinstance(X, pd.DataFrame):
        raise ValueError("X must be DataFrames")

    if not isinstance(F, pd.DataFrame) and not isinstance(B, pd.DataFrame):
        raise ValueError("F and B must be DataFrames")

    assets = X.columns.tolist()

    if B is not None:
        B = np.array(B, ndmin=2)
        if const == True:
            alpha = B[:, :1]
            B = B[:, 1:]

    mu_f = np.array(mean_vector(F, method=method_mu, **dict_mu), ndmin=2)
    mu_f = (mu_f - rf).T

    tau = 1 / X.shape[0]

    S_f = np.array(covar_matrix(F, method=method_cov, **dict_cov), ndmin=2)
    S = B @ S_f @ B.T

    D = X.to_numpy() - F @ B.T
    D = np.diag(D.var())
    S = S + D

    Omega_f = np.array(np.diag(np.diag(P_f @ (tau * S_f) @ P_f.T)), ndmin=2)

    S_hat = inv(inv(S_f) + P_f.T @ inv(Omega_f) @ P_f)

    Pi_hat = S_hat @ (inv(S_f) @ mu_f + P_f.T @ inv(Omega_f) @ Q_f)

    S_blb = inv(inv(S) - inv(S) @ B @ inv(inv(S_hat) + B.T @ inv(S) @ B) @ B.T @ inv(S))

    Pi_blb = (
        S_blb @ inv(S) @ B @ inv(inv(S_hat) + B.T @ inv(S) @ B) @ inv(S_hat) @ Pi_hat
    )

    mu = Pi_blb + rf

    if const == True:
        mu = mu + alpha
    mu = mu.T
    cov = S_blb
    w = inv(delta * cov) @ mu.T

    mu = pd.DataFrame(mu, columns=assets)
    cov = pd.DataFrame(cov, index=assets, columns=assets)
    w = pd.DataFrame(w, index=assets)

    return mu, cov, w


def bootstrapping(
    X,
    kind="stationary",
    q=0.05,
    n_sim=6000,
    window=3,
    diag=False,
    threshold=1e-15,
    seed=0,
):
    r"""
    Estimates the uncertainty sets of mean and covariance matrix through the selected
    bootstrapping method.

    Parameters
    ----------
    X : DataFrame of shape (n_samples, n_assets)
        Assets returns DataFrame, where n_samples is the number of
        observations and n_assets is the number of assets.
    kind : str
        The bootstrapping method. The default value is 'stationary'. Possible values are:

        - 'stationary': stationary bootstrapping method, see `StationaryBootstrap <https://bashtage.github.io/arch/bootstrap/generated/arch.bootstrap.StationaryBootstrap.html#arch.bootstrap.StationaryBootstrap>`_ for more details.
        - 'circular': circular bootstrapping method, see `CircularBlockBootstrap <https://bashtage.github.io/arch/bootstrap/generated/arch.bootstrap.CircularBlockBootstrap.html#arch.bootstrap.CircularBlockBootstrap>`_ for more details.
        - 'moving': moving bootstrapping method, see `MovingBlockBootstrap <https://bashtage.github.io/arch/bootstrap/generated/arch.bootstrap.MovingBlockBootstrap.html#arch.bootstrap.MovingBlockBootstrap>`_ for more details.
    q : scalar
        Significance level for box and elliptical constraints.
        The default is 0.05.
    n_sim : scalar
        Number of simulations of the bootstrapping method.
        The default is 6000.
    window: int
        Block size of the bootstrapping method. Must be greather than 1
        and lower than the n_samples - n_factors + 1
        The default is 3.
    diag: bool
        If consider only the main diagonal of covariance matrices of estimation
        errors following :cite:`b-fabozzi2007robust`. The default is False.
    threshold: float
        Parameter used to fix covariance matrices in case they are not positive semidefinite.
        The default is 1e-15.
    seed: int
        Seed used to generate random numbers for bootstrapping method.
        The default is 0.

    Returns
    -------
    mu_l : DataFrame
        The q/2 percentile of mean vector obtained through the selected
        bootstrapping method.
    mu_u : DataFrame
        The 1-q/2 percentile of mean vector obtained through the selected
        bootstrapping method.
    cov_l : DataFrame
        The q/2 percentile of covariance matrix obtained through the selected
        bootstrapping method.
    cov_u : DataFrame
        The 1-q/2 percentile of covariance matrix obtained through the selected
        bootstrapping method.
    cov_mu : DataFrame
        The covariance matrix of estimation errors of mean vector obtained
        through the selected bootstrapping method.
    cov_sigma : DataFrame
        The covariance matrix of estimation errors of covariance matrix
        obtained through the selected bootstrapping method.
    k_mu : DataFrame
        The square root of size of elliptical constraint of mean vector
        estimation error based on 1-q percentile.
    k_sigma : DataFrame
        The square root of size of elliptical constraint of covariance matrix
        estimation error based on 1-q percentile.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    """

    if not isinstance(X, pd.DataFrame):
        raise ValueError("X must be a DataFrame")

    if window >= X.shape[0] - window + 1:
        raise ValueError("block must be lower than  n_samples - window + 1")
    elif window <= 1:
        raise ValueError("block must be greather than 1")

    cols = X.columns.tolist()
    cols_2 = [i + "-" + j for i in cols for j in cols]
    T, n = X.shape

    mu = X.mean().to_numpy().reshape(1, n)
    vec_Sigma = X.cov().to_numpy().reshape((1, n**2), order="F")

    mus = np.zeros((n_sim, 1, n))
    covs = np.zeros((n_sim, n, n))

    if kind == "stationary":
        gen = bs.StationaryBootstrap(window, X, seed=seed)
    elif kind == "circular":
        gen = bs.CircularBlockBootstrap(window, X, seed=seed)
    elif kind == "moving":
        gen = bs.MovingBlockBootstrap(window, X, seed=seed)
    else:
        raise ValueError("kind only can be 'stationary', 'circular' or 'moving'")

    i = 0
    for data in gen.bootstrap(n_sim):
        A = data[0][0]
        mus[i] = A.mean().to_numpy().reshape(1, n)
        covs[i] = A.cov().to_numpy()
        i += 1

    # Box Constraint for Mean
    mu_l = np.percentile(mus, q=q / 2 * 100, axis=0, keepdims=True).reshape(1, n)
    mu_u = np.percentile(mus, q=(1 - q / 2) * 100, axis=0, keepdims=True).reshape(1, n)
    mu_l = pd.DataFrame(mu_l, index=[0], columns=cols)
    mu_u = pd.DataFrame(mu_u, index=[0], columns=cols)

    # Box Constraint for Covariance
    cov_l = np.percentile(covs, q=q / 2 * 100, axis=0, keepdims=True).reshape(n, n)
    cov_u = np.percentile(covs, q=(1 - q / 2) * 100, axis=0, keepdims=True).reshape(
        n, n
    )
    cov_l = pd.DataFrame(cov_l, index=cols, columns=cols)
    cov_u = pd.DataFrame(cov_u, index=cols, columns=cols)

    # Check and fix if upper and lower bound for Covariance are positive
    # semidefinite and fix when they are not
    if af.is_pos_def(cov_l) == False:
        cov_l = af.cov_fix(cov_l, method="clipped", threshold=threshold)
    if af.is_pos_def(cov_u) == False:
        cov_u = af.cov_fix(cov_u, method="clipped", threshold=threshold)

    # Elliptical Constraint for Mean
    A_mu = mus.reshape(n_sim, n) - np.repeat(mu, n_sim, axis=0)
    cov_mu = np.cov(A_mu, rowvar=False)
    if diag == True:
        cov_mu = np.diag(np.diag(cov_mu))
    k_mus = np.diag(A_mu @ inv(cov_mu) @ A_mu.T)
    k_mu = np.percentile(k_mus, q=(1 - q) * 100) ** 0.5
    cov_mu = pd.DataFrame(cov_mu, index=cols, columns=cols)

    # Elliptical Constraint for Covariance
    A_Sigma = covs.reshape((n_sim, n**2), order="F")
    A_Sigma = A_Sigma - np.repeat(vec_Sigma, n_sim, axis=0)
    cov_sigma = np.cov(A_Sigma, rowvar=False)
    cov_sigma = af.cov_fix(cov_sigma, method="clipped", threshold=threshold)
    if diag == True:
        cov_sigma = np.diag(np.diag(cov_sigma))
    if af.is_pos_def(cov_sigma) == False:
        cov_sigma = af.cov_fix(cov_sigma, method="clipped", threshold=threshold)
    k_sigmas = np.diag(A_Sigma @ inv(cov_sigma) @ A_Sigma.T)
    k_sigma = np.percentile(k_sigmas, q=(1 - q) * 100) ** 0.5
    cov_sigma = pd.DataFrame(cov_sigma, index=cols_2, columns=cols_2)

    return mu_l, mu_u, cov_l, cov_u, cov_mu, cov_sigma, k_mu, k_sigma


def normal_simulation(X, q=0.05, n_sim=6000, diag=False, threshold=1e-15, seed=0):
    r"""
    Estimates the uncertainty sets of mean and covariance matrix assuming that
    assets returns follows a multivariate normal distribution.

    Parameters
    ----------
    X : DataFrame of shape (n_samples, n_assets)
        Assets returns DataFrame, where n_samples is the number of
        observations and n_assets is the number of assets.
    q : scalar
        Significance level for box and elliptical constraints.
        The default is 0.05.
    n_sim : scalar
        Number of simulations of the bootstrapping method.
        The default is 6000.
    diag: bool
        If consider only the main diagonal of covariance matrices of estimation
        errors following :cite:`b-fabozzi2007robust`. The default is False.
    threshold: float
        Parameter used to fix covariance matrices in case they are not positive
        semidefinite. The default is 1e-10.
    seed: int
        Seed used to generate random numbers for simulation.
        The default is 0.

    Returns
    -------
    mu_l : DataFrame
        The q/2 percentile of mean vector obtained through the normal
        simulation.
    mu_u : DataFrame
        The 1-q/2 percentile of mean vector obtained through the normal
        simulation.
    cov_l : DataFrame
        The q/2 percentile of covariance matrix obtained through the normal
        simulation.
    cov_u : DataFrame
        The 1-q/2 percentile of covariance matrix obtained through the normal
        simulation.
    cov_mu : DataFrame
        The covariance matrix of estimation errors of mean vector obtained
        through the normal simulation.
    cov_sigma : DataFrame
        The covariance matrix of estimation errors of covariance matrix
        obtained through the normal simulation.
    k_mu : DataFrame
        The square root of size of elliptical constraint of mean vector
        estimation error based on 1-q percentile.
    k_sigma : DataFrame
        The square root of size of elliptical constraint of covariance matrix
        estimation error based on 1-q percentile.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    """

    if not isinstance(X, pd.DataFrame):
        raise ValueError("X must be a DataFrame")

    cols = X.columns.tolist()
    cols_2 = [i + "-" + j for i in cols for j in cols]
    T, n = X.shape

    # Set initial parameters based on assumption of normality
    mu = X.mean().to_numpy().reshape(1, n)
    vec_Sigma = X.cov().to_numpy().reshape((1, n**2), order="F")
    Sigma = X.cov().to_numpy()
    cov_mu = Sigma / T
    K = cf.commutation_matrix(T=n, n=n)
    I = np.identity(n**2)
    cov_sigma = T * (I + K) @ np.kron(cov_mu, cov_mu)
    if diag == True:
        cov_sigma = np.diag(np.diag(cov_sigma))
    if af.is_pos_def(cov_sigma) == False:
        cov_sigma = af.cov_fix(cov_sigma, method="clipped", threshold=threshold)
    cov_sigma = pd.DataFrame(cov_sigma, index=cols_2, columns=cols_2)

    # Box Constraint for Mean
    delta_mu = st.norm.ppf(1 - q / 2) * np.sqrt(np.diag(cov_mu)).reshape(-1, 1)
    mu_l = mu - delta_mu.T
    mu_u = mu + delta_mu.T
    mu_l = pd.DataFrame(mu_l, index=[0], columns=cols)
    mu_u = pd.DataFrame(mu_u, index=[0], columns=cols)

    # Box Constraints for Covariance
    rs = np.random.RandomState(seed=seed)
    covs = st.wishart.rvs(T, cov_mu, size=n_sim, random_state=rs)
    cov_l = np.percentile(covs, q=q / 2, axis=0)
    cov_u = np.percentile(covs, q=1 - q / 2, axis=0)
    cov_l = pd.DataFrame(cov_l, index=cols, columns=cols)
    cov_u = pd.DataFrame(cov_u, index=cols, columns=cols)

    # Check and fix if upper and lower bound for Covariance are positive
    # semidefinite and fix when they are not
    if af.is_pos_def(cov_l) == False:
        cov_l = af.cov_fix(cov_l, method="clipped", threshold=threshold)
    if af.is_pos_def(cov_u) == False:
        cov_u = af.cov_fix(cov_u, method="clipped", threshold=threshold)

    # Elliptical Constraint for Mean
    A_mu = rs.multivariate_normal(mu.ravel(), cov_mu, size=n_sim)
    # cov_mu =  np.cov(A_mu - np.repeat(mu, n_sim, axis=0), rowvar=False)
    if diag == True:
        cov_mu = np.diag(np.diag(cov_mu))
    k_mus = np.diag(A_mu @ inv(cov_mu) @ A_mu.T)
    k_mu = np.percentile(k_mus, q=1 - q) ** 0.5
    # k_mu = st.chi2.ppf(1 - q, df=n) ** 0.5
    cov_mu = pd.DataFrame(cov_mu, index=cols, columns=cols)

    # Elliptical Constraint for Covariance
    A_Sigma = covs.reshape((n_sim, n**2), order="F")
    A_Sigma = A_Sigma - np.repeat(vec_Sigma, n_sim, axis=0)
    A_cov_sigma = np.cov(A_Sigma, rowvar=False)
    if diag == True:
        A_cov_sigma = np.diag(np.diag(A_cov_sigma))
    if af.is_pos_def(A_cov_sigma) == False:
        A_cov_sigma = af.cov_fix(A_cov_sigma, method="clipped", threshold=threshold)
    k_sigmas = np.diag(A_Sigma @ inv(A_cov_sigma) @ A_Sigma.T)
    k_sigma = np.percentile(k_sigmas, q=1 - q) ** 0.5

    return mu_l, mu_u, cov_l, cov_u, cov_mu, cov_sigma, k_mu, k_sigma
