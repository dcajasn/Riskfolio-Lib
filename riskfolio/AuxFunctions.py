""""""  #
"""
Copyright (c) 2020-2021, Dany Cajas
All rights reserved.
This work is licensed under BSD 3-Clause "New" or "Revised" License.
License available at https://github.com/dcajasn/Riskfolio-Lib/blob/master/LICENSE.txt
"""

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import linalg as LA
from statsmodels.stats.correlation_tools import cov_nearest
from scipy.sparse import csr_matrix
from scipy.spatial.distance import pdist, squareform
import scipy.cluster.hierarchy as hr
import scipy.stats as st
from sklearn.metrics import mutual_info_score
from astropy.stats import knuth_bin_width, freedman_bin_width, scott_bin_width

###############################################################################
# Aditional Matrix Functions
###############################################################################


def is_pos_def(cov, threshold=1e-8):
    r"""
    Indicate if a matrix is positive (semi)definite.

    Parameters
    ----------
    cov : ndarray
        Features covariance matrix of shape (n_features, n_features), where
        n_features is the number of features.

    Returns
    -------
    value : bool
        True if matrix is positive (semi)definite.

    Raises
    ------
        ValueError when the value cannot be calculated.

    """
    cov_ = np.array(cov, ndmin=2)
    w, V = LA.eigh(cov_, lower=True, check_finite=True)
    value = np.all(w >= threshold)

    return value


def correl_matrix(cov):
    r"""
    Generate a correlation matrix from a covariance matrix cov.

    Parameters
    ----------
    cov : ndarray
        Assets covariance matrix of shape n_features x n_features, where
        n_features is the number of features.

    Returns
    -------
    corr : ndarray
        A correlation matrix.

    Raises
    ------
        ValueError when the value cannot be calculated.

    """
    flag = False
    if isinstance(cov, pd.DataFrame):
        cols = cov.columns.tolist()
        flag = True

    cov1 = np.array(cov, ndmin=2)
    corr = np.array(cov, ndmin=2)
    m, n = cov.shape
    for i in range(0, m):
        for j in range(0, n):
            corr[i, j] = cov1[i, j] / np.sqrt(cov1[i, i] * cov1[j, j])

    if flag:
        corr = pd.DataFrame(corr, index=cols, columns=cols)

    return corr


def cov_fix(cov, method="clipped", **kwargs):
    r"""
    Fix a covariance matrix to a positive definite matrix.

    Parameters
    ----------
    cov : ndarray
        Features covariance matrix of shape n_features x n_features, where
        n_features is the number of features.
    method : str
        The default value is 'clipped', see more in `cov_nearest <https://www.statsmodels.org/stable/generated/statsmodels.stats.correlation_tools.cov_nearest.html>`_.
    **kwargs
        Other parameters from `cov_nearest <https://www.statsmodels.org/stable/generated/statsmodels.stats.correlation_tools.cov_nearest.html>`_.

    Returns
    -------
    cov_ : bool
        A positive definite covariance matrix.

    Raises
    ------
        ValueError when the value cannot be calculated.

    """
    flag = False
    if isinstance(cov, pd.DataFrame):
        cols = cov.columns.tolist()
        flag = True

    cov_ = np.array(cov, ndmin=2)
    cov_ = cov_nearest(cov_, method=method, **kwargs)
    cov_ = np.array(cov_, ndmin=2)

    if flag:
        cov_ = pd.DataFrame(cov_, index=cols, columns=cols)

    return cov_


def cov_returns(cov, seed=0):
    r"""
    Generate a matrix of returns that have a covariance matrix cov.

    Parameters
    ----------
    cov : ndarray
        Assets covariance matrix of shape n_features x n_features, where
        n_features is the number of features.

    Returns
    -------
    a : ndarray
        A matrix of returns that have a covariance matrix cov.

    Raises
    ------
        ValueError when the value cannot be calculated.

    """

    rs = np.random.RandomState(seed)
    n = len(cov)
    a = np.array(rs.randn(n + 10, n), ndmin=2)

    for i in range(0, 5):
        cov_ = np.cov(a.T)
        L = np.array(np.linalg.cholesky(cov_), ndmin=2)
        a = a @ np.linalg.inv(L).T
        cov_ = np.cov(a.T)
        desv_ = np.sqrt(np.array(np.diag(cov_), ndmin=2))
        a = (np.array(a) - np.mean(a, axis=0)) / np.array(desv_)

    L1 = np.array(np.linalg.cholesky(cov), ndmin=2)
    a = a @ L1.T

    return a


def commutation_matrix(cov):
    r"""
    Generate the commutation matrix of the covariance matrix cov.

    Parameters
    ----------
    cov : ndarray
        Assets covariance matrix of shape n_features x n_features, where
        n_features is the number of features.

    Returns
    -------
    K : ndarray
        The commutation matrix of the covariance matrix cov.

    Raises
    ------
        ValueError when the value cannot be calculated.

    """
    (m, n) = cov.shape
    row = np.arange(m * n)
    col = row.reshape((m, n), order="F").ravel()
    data = np.ones(m * n, dtype=np.int8)
    K = csr_matrix((data, (row, col)), shape=(m * n, m * n))
    K = K.toarray()

    return K


###############################################################################
# Aditional Codependence Functions
###############################################################################


def dcorr(X, Y):
    r"""
    Calculate the distance correlation between two variables :cite:`d-Szekely`.

    Parameters
    ----------
    X : 1d-array
        Returns series, must have of shape n_sample x 1.
    Y : 1d-array
        Returns series, must have of shape n_sample x 1.

    Returns
    -------
    value : float
        The distance correlation between variables X and Y.

    Raises
    ------
        ValueError when the value cannot be calculated.

    """

    X = np.atleast_1d(X)
    Y = np.atleast_1d(Y)

    if np.prod(X.shape) == len(X):
        X = X[:, None]
    if np.prod(Y.shape) == len(Y):
        Y = Y[:, None]

    X = np.atleast_2d(X)
    Y = np.atleast_2d(Y)
    n = X.shape[0]

    if Y.shape[0] != X.shape[0]:
        raise ValueError("Number of samples must match")

    a = squareform(pdist(X))
    b = squareform(pdist(Y))
    A = a - a.mean(axis=0)[None, :] - a.mean(axis=1)[:, None] + a.mean()
    B = b - b.mean(axis=0)[None, :] - b.mean(axis=1)[:, None] + b.mean()

    dcov2_xy = (A * B).sum() / float(n * n)
    dcov2_xx = (A * A).sum() / float(n * n)
    dcov2_yy = (B * B).sum() / float(n * n)
    value = np.sqrt(dcov2_xy) / np.sqrt(np.sqrt(dcov2_xx) * np.sqrt(dcov2_yy))

    return value


def dcorr_matrix(X):
    r"""
    Calculate the distance correlation matrix of n variables.

    Parameters
    ----------
    X : ndarray or
        Returns series of shape n_sample x n_features.

    Returns
    -------
    corr : ndarray
        The distance correlation matrix of shape n_features x n_features.

    Raises
    ------
        ValueError when the value cannot be calculated.

    """
    flag = False
    if isinstance(X, pd.DataFrame):
        cols = X.columns.tolist()
        X1 = X.to_numpy()
        flag = True
    else:
        X1 = X.copy()

    n = X1.shape[1]
    corr = np.ones((n, n))
    indices = np.triu_indices(n, 1)

    for i, j in zip(indices[0], indices[1]):
        corr[i, j] = dcorr(X1[:, i], X1[:, j])
        corr[j, i] = corr[i, j]

    if flag:
        corr = pd.DataFrame(corr, index=cols, columns=cols)
    else:
        corr = pd.DataFrame(corr)

    return corr


def numBins(n_samples, corr=None):
    r"""
    Calculate the optimal number of bins for discretization of mutual
    information and variation of information.

    Parameters
    ----------
    n_samples : integer
        Number of samples.

    corr : float, optional
        Correlation coefficient of variables. The default value is None.

    Returns
    -------
    bins : int
        The optimal number of bins.

    Raises
    ------
        ValueError when the value cannot be calculated.

    """
    # univariate case
    if corr is None:
        z = (
            8 + 324 * n_samples + 12 * (36 * n_samples + 729 * n_samples ** 2) ** 0.5
        ) ** (1 / 3)
        b = np.round(z / 6 + 2 / (3 * z) + 1 / 3)
    # bivariate case
    else:
        b = np.round(
            2 ** -0.5 * (1 + (1 + 24 * n_samples / (1 - corr ** 2)) ** 0.5) ** 0.5
        )

    bins = np.int32(b)

    return bins


def mutual_info_matrix(X, bins_info="KN", normalize=True):
    r"""
    Calculate the mutual information matrix of n variables.

    Parameters
    ----------
    X : ndarray
        Returns series of shape n_sample x n_features.
    bins_info: int or str
        Number of bins used to calculate mutual information. The default
        value is 'KN'. Posible values are:

        - 'KN': Knuth's choice method. See more in `knuth_bin_width <https://docs.astropy.org/en/stable/api/astropy.stats.knuth_bin_width.html>`_.
        - 'FD': Freedman–Diaconis' choice method. See more in `freedman_bin_width <https://docs.astropy.org/en/stable/api/astropy.stats.freedman_bin_width.html>`_.
        - 'SC': Scotts' choice method. See more in `scott_bin_width <https://docs.astropy.org/en/stable/api/astropy.stats.scott_bin_width.html>`_.
        - 'HGR': Hacine-Gharbi and Ravier' choice method.
        - int: integer value choice by user.

    normalize: bool
        If normalize variation of information. The default value is True.

    Returns
    -------
    corr : ndarray
        The mutual information matrix of shape n_features x n_features.

    Raises
    ------
        ValueError when the value cannot be calculated.

    """
    flag = False
    if isinstance(X, pd.DataFrame):
        cols = X.columns.tolist()
        X1 = X.to_numpy()
        flag = True
    else:
        X1 = X.copy()

    m = X1.shape[0]
    n = X1.shape[1]
    mat = np.zeros((n, n))
    indices = np.triu_indices(n)

    for i, j in zip(indices[0], indices[1]):
        if bins_info == "KN":
            k1 = (np.max(X1[:, i]) - np.min(X1[:, i])) / knuth_bin_width(X1[:, i])
            bins = np.int32(np.round(k1))
            if i != j:
                k2 = (np.max(X1[:, j]) - np.min(X1[:, j])) / knuth_bin_width(X1[:, j])
                bins = np.int32(np.round(np.maximum(k1, k2)))
        elif bins_info == "FD":
            k1 = (np.max(X1[:, i]) - np.min(X1[:, i])) / freedman_bin_width(X1[:, i])
            bins = np.int32(np.round(k1))
            if i != j:
                k2 = (np.max(X1[:, j]) - np.min(X1[:, j])) / freedman_bin_width(
                    X1[:, j]
                )
                bins = np.int32(np.round(np.maximum(k1, k2)))
        elif bins_info == "SC":
            k1 = (np.max(X1[:, i]) - np.min(X1[:, i])) / scott_bin_width(X1[:, i])
            bins = np.int32(np.round(k1))
            if i != j:
                k2 = (np.max(X1[:, j]) - np.min(X1[:, j])) / scott_bin_width(X1[:, j])
                bins = np.int32(np.round(np.maximum(k1, k2)))
        elif bins_info == "HGR":
            corr = np.corrcoef(X1[:, i], X1[:, j])[0, 1]
            if corr == 1:
                bins = numBins(m, None)
            else:
                bins = numBins(m, corr)
        elif isinstance(bins_info, np.int32) or isinstance(bins_info, int):
            bins = bins_info

        cXY = np.histogram2d(X1[:, i], X1[:, j], bins)[0]
        hX = st.entropy(np.histogram(X1[:, i], bins)[0])  # marginal
        hY = st.entropy(np.histogram(X1[:, j], bins)[0])  # marginal
        iXY = mutual_info_score(None, None, contingency=cXY)  # mutual information
        if normalize == True:
            iXY = iXY / np.min([hX, hY])  # normalized mutual information
            # hXY = hX + hY - iXY # joint
            # hX_Y = hXY - hY # conditional
            # hY_X = hXY - hX # conditional

        mat[i, j] = iXY
        mat[j, i] = mat[i, j]

    mat = np.clip(np.round(mat, 8), a_min=0.0, a_max=np.inf)

    if flag:
        mat = pd.DataFrame(mat, index=cols, columns=cols)
    else:
        mat = pd.DataFrame(mat)

    return mat


def var_info_matrix(X, bins_info="KN", normalize=True):
    r"""
    Calculate the variation of information matrix of n variables.

    Parameters
    ----------
    X : ndarray
        Returns series of shape n_sample x n_features.
    bins_info: int or str
        Number of bins used to calculate variation of information. The default
        value is 'KN'. Posible values are:

        - 'KN': Knuth's choice method. See more in `knuth_bin_width <https://docs.astropy.org/en/stable/api/astropy.stats.knuth_bin_width.html>`_.
        - 'FD': Freedman–Diaconis' choice method. See more in `freedman_bin_width <https://docs.astropy.org/en/stable/api/astropy.stats.freedman_bin_width.html>`_.
        - 'SC': Scotts' choice method. See more in `scott_bin_width <https://docs.astropy.org/en/stable/api/astropy.stats.scott_bin_width.html>`_.
        - 'HGR': Hacine-Gharbi and Ravier' choice method.
        - int: integer value choice by user.

    normalize: bool
        If normalize variation of information. The default value is True.

    Returns
    -------
    corr : ndarray
        The mutual information matrix of shape n_features x n_features.

    Raises
    ------
        ValueError when the value cannot be calculated.

    """
    flag = False
    if isinstance(X, pd.DataFrame):
        cols = X.columns.tolist()
        X1 = X.to_numpy()
        flag = True
    else:
        X1 = X.copy()

    m = X1.shape[0]
    n = X1.shape[1]
    mat = np.zeros((n, n))
    indices = np.triu_indices(n)

    for i, j in zip(indices[0], indices[1]):
        if bins_info == "KN":
            k1 = (np.max(X1[:, i]) - np.min(X1[:, i])) / knuth_bin_width(X1[:, i])
            bins = np.int32(np.round(k1))
            if i != j:
                k2 = (np.max(X1[:, j]) - np.min(X1[:, j])) / knuth_bin_width(X1[:, j])
                bins = np.int32(np.round(np.maximum(k1, k2)))
        elif bins_info == "FD":
            k1 = (np.max(X1[:, i]) - np.min(X1[:, i])) / freedman_bin_width(X1[:, i])
            bins = np.int32(np.round(k1))
            if i != j:
                k2 = (np.max(X1[:, j]) - np.min(X1[:, j])) / freedman_bin_width(
                    X1[:, j]
                )
                bins = np.int32(np.round(np.maximum(k1, k2)))
        elif bins_info == "SC":
            k1 = (np.max(X1[:, i]) - np.min(X1[:, i])) / scott_bin_width(X1[:, i])
            bins = np.int32(np.round(k1))
            if i != j:
                k2 = (np.max(X1[:, j]) - np.min(X1[:, j])) / scott_bin_width(X1[:, j])
                bins = np.int32(np.round(np.maximum(k1, k2)))
        elif bins_info == "HGR":
            corr = np.corrcoef(X1[:, i], X1[:, j])[0, 1]
            if corr == 1:
                bins = numBins(m, None)
            else:
                bins = numBins(m, corr)
        elif isinstance(bins_info, np.int32) or isinstance(bins_info, int):
            bins = bins_info

        cXY = np.histogram2d(X1[:, i], X1[:, j], bins)[0]
        hX = st.entropy(np.histogram(X1[:, i], bins)[0])  # marginal
        hY = st.entropy(np.histogram(X1[:, j], bins)[0])  # marginal
        iXY = mutual_info_score(None, None, contingency=cXY)  # mutual information
        vXY = hX + hY - 2 * iXY  # variation of information
        if normalize == True:
            hXY = hX + hY - iXY  # joint
            vXY = vXY / hXY  # normalized variation of information

        mat[i, j] = vXY
        mat[j, i] = mat[i, j]

    mat = np.clip(np.round(mat, 8), a_min=0.0, a_max=np.inf)

    if flag:
        mat = pd.DataFrame(mat, index=cols, columns=cols)
    else:
        mat = pd.DataFrame(mat)

    return mat


def ltdi_matrix(X, alpha=0.05):
    r"""
    Calculate the lower tail dependence index matrix using the empirical
    approach.

    Parameters
    ----------
    X : ndarray
        Returns series of shape n_sample x n_features.
    alpha : float, optional
        Significance level for lower tail dependence index.
        The default is 0.05.

    Returns
    -------
    corr : ndarray
        The lower tail dependence index matrix of shape n_features x
        n_features.

    Raises
    ------
        ValueError when the value cannot be calculated.

    """

    flag = False
    if isinstance(X, pd.DataFrame):
        cols = X.columns.tolist()
        X1 = X.to_numpy()
        flag = True
    else:
        X1 = X.copy()

    m = X1.shape[0]
    n = X1.shape[1]
    k = np.int(np.ceil(m * alpha))
    mat = np.ones((n, n))

    if k > 0:
        indices = np.triu_indices(n)

        for i, j in zip(indices[0], indices[1]):
            u = np.sort(X1[:, i])[k - 1]
            v = np.sort(X1[:, j])[k - 1]
            ltd = (
                np.sum(np.where(np.logical_and(X1[:, i] <= u, X1[:, j] <= v), 1, 0)) / k
            )

            mat[i, j] = ltd
            mat[j, i] = mat[i, j]

        for i in range(0, n):
            u = np.sort(X1[:, i])[k - 1]
            v = np.sort(X1[:, i])[k - 1]
            ltd = (
                np.sum(np.where(np.logical_and(X1[:, i] <= u, X1[:, i] <= v), 1, 0)) / k
            )

            mat[i, i] = ltd

    mat = np.clip(np.round(mat, 8), a_min=1.0e-8, a_max=1)
    
    if flag:
        mat = pd.DataFrame(mat, index=cols, columns=cols)
    else:
        mat = pd.DataFrame(mat)

    return mat


def two_diff_gap_stat(codep, dist, clusters, max_k=10):
    r"""
    Calculate the optimal number of clusters based on the two difference gap
    statistic :cite:`d-twogap`.

    Parameters
    ----------
    codep : DataFrame
        A codependence matrix.
    dist : str, optional
        A distance measure based on the codependence matrix.
    clusters : string, optional
        The hierarchical clustering encoded as a linkage matrix, see `linkage <https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html?highlight=linkage#scipy.cluster.hierarchy.linkage>`_ for more details.
    max_k : int, optional
        Max number of clusters used by the two difference gap statistic
        to find the optimal number of clusters. The default is 10.

    Returns
    -------
    k : int
        The optimal number of clusters based on the two difference gap statistic.

    Raises
    ------
        ValueError when the value cannot be calculated.

    """
    # cluster levels over from 1 to N-1 clusters
    cluster_lvls = pd.DataFrame(hr.cut_tree(clusters), index=codep.columns)
    num_k = cluster_lvls.columns  # save column with number of clusters
    cluster_lvls = cluster_lvls.iloc[:, ::-1]  # reverse order to start with 1 cluster
    cluster_lvls.columns = num_k  # set columns to number of cluster
    W_list = []

    # get within-cluster dissimilarity for each k
    for k in range(min(len(cluster_lvls.columns), max_k)):
        level = cluster_lvls.iloc[:, k]  # get k clusters
        D_list = []  # within-cluster distance list

        for i in range(np.max(level.unique()) + 1):
            cluster = level.loc[level == i]
            # Based on correlation distance
            cluster_dist = dist.loc[cluster.index, cluster.index]  # get distance
            cluster_pdist = squareform(cluster_dist, checks=False)
            if cluster_pdist.shape[0] != 0:
                D = np.nan_to_num(cluster_pdist.mean())
                D_list.append(D)  # append to list

        W_k = np.sum(D_list)
        W_list.append(W_k)

    W_list = pd.Series(W_list)
    n = codep.shape[0]
    limit_k = int(min(max_k, np.sqrt(n)))
    gaps = W_list.shift(2) + W_list - 2 * W_list.shift(1)
    gaps = gaps[0:limit_k]
    if gaps.isna().all():
        k = len(gaps)
    else:
        k = int(gaps.idxmax() + 2)

    return k


###############################################################################
# Other Aditional Functions
###############################################################################


def round_values(data, decimals=4, wider=False):
    """
    This function help us to round values to values close or away from zero.

    Parameters
    ----------
    data : np.ndarray, pd.Series or pd.DataFrame
        Data that are going to be rounded.
    decimals : integer
        Number of decimals to round.
    wider : float
        False if round to values close to zero, True if round to values away
        from zero.

    Returns
    -------
    value : np.ndarray, pd.Series or pd.DataFrame
        Data rounded using selected method.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    """

    if wider == True:
        value = np.where(
            data >= 0,
            np.ceil(data * 10 ** decimals) / 10 ** decimals,
            np.floor(data * 10 ** decimals) / 10 ** decimals,
        )
    elif wider == False:
        value = np.where(
            data >= 0,
            np.floor(data * 10 ** decimals) / 10 ** decimals,
            np.ceil(data * 10 ** decimals) / 10 ** decimals,
        )

    if isinstance(data, pd.DataFrame):
        value = pd.DataFrame(value, columns=data.columns, index=data.index)
    if isinstance(data, pd.Series):
        value = pd.Series(value, index=data.index)

    return value


def weights_discretizetion(
    weights, prices, capital=1000000, w_decimal=6, ascending=False
):
    """
    This function help us to find the number of shares that must be bought or
    sold to achieve portfolio weights according the prices of assets and the
    invested capital.

    Parameters
    ----------
    weights : pd.Series or pd.DataFrame
        Vector of weights of size n_assets x 1.
    prices : pd.Series or pd.DataFrame
        Vector of prices of size n_assets x 1.
    capital : float, optional
        Capital invested. The default value is 1000000.
    w_decimal : int, optional
        Number of decimals use to round the portfolio weights. The default
        value is 6.
    ascending : bool, optional
        If True assigns excess capital to assets with lower weights, else,
        to assets with higher weights. The default value is False.

    Returns
    -------
    n_shares : pd.DataFrame
        Number of shares that must be bought or sold to achieve portfolio
        weights.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    """

    if isinstance(weights, pd.Series):
        w = weights.to_frame().copy()
    elif isinstance(weights, pd.DataFrame):
        if weights.shape[0] == 1:
            w = weights.T.copy()
        elif weights.shape[1] == 1:
            w = weights.copy()
            pass
        else:
            raise ValueError("weights must have size n_assets x 1")
    else:
        raise ValueError("weights must be DataFrame")

    if isinstance(prices, pd.Series):
        p = prices.to_frame().copy()
    elif isinstance(prices, pd.DataFrame):
        if prices.shape[0] == 1:
            p = prices.T.copy()
        elif prices.shape[1] == 1:
            p = prices.copy()
            pass
        else:
            raise ValueError("prices must have size n_assets x 1")
    else:
        raise ValueError("prices must be DataFrame")

    w.columns = [0]
    p.columns = [0]

    total = w.sum().item()
    w = round_values(w, decimals=w_decimal, wider=False)
    w.loc[w.idxmin().tolist()] = w.loc[w.idxmin().tolist()] + (total - w.sum()).item()

    n_shares = round_values(capital * w / p, decimals=0, wider=False)

    excedent = [capital + 1, capital]
    i = 1
    while excedent[i] < excedent[i - 1]:
        new_capital = (n_shares.T @ p).iloc[0, 0]
        excedent.append(capital - new_capital)
        new_shares = round_values(excedent[-1] * w / p, 0)
        n_shares += new_shares
        i += 1

    n_shares_1 = capital * w / p

    excedent = capital - (n_shares.T @ p).iloc[0, 0]
    i = 1

    d_shares = np.abs(n_shares_1) - np.abs(n_shares)
    d_shares = np.where(d_shares > 0, n_shares_1 - n_shares, 0)
    d_shares = round_values(d_shares, decimals=0, wider=True)
    d_shares = pd.DataFrame(d_shares, columns=w.columns, index=w.index)

    n_shares_1 = capital * w / p

    excedent = capital - (n_shares.T @ p).iloc[0, 0]

    d_shares = np.abs(n_shares_1) - np.abs(n_shares)
    d_shares = np.where(d_shares > 0, n_shares_1 - n_shares, 0)
    d_shares = round_values(d_shares, decimals=0, wider=True)
    d_shares = pd.DataFrame(d_shares, columns=w.columns, index=w.index)

    order = w.sort_values(by=0, ascending=ascending).index.tolist()
    d_list = d_shares[d_shares[0] == 1].index.tolist()

    for i in order:
        if i in d_list:
            new_shares = round_values(excedent / p.loc[i, 0], 0).item()
            if new_shares > 0:
                n_shares.loc[i] += new_shares
                excedent = capital - (n_shares.T @ p).iloc[0, 0]

    return n_shares


def color_list(k):
    r"""
    This function creates a list of colors.

    Parameters
    ----------
    k : int
        Number of colors.

    Returns
    -------
    colors : list
        A list of colors.
    """

    colors = []

    if k <= 10:
        for i in range(10):
            colors.append(mpl.colors.rgb2hex(plt.get_cmap("tab10").colors[i]))
    elif k <= 20:
        for i in range(20):
            colors.append(mpl.colors.rgb2hex(plt.get_cmap("tab20").colors[i]))
    elif k <= 40:
        for i in range(20):
            colors.append(mpl.colors.rgb2hex(plt.get_cmap("tab20").colors[i]))
        for i in range(20):
            colors.append(mpl.colors.rgb2hex(plt.get_cmap("tab20b").colors[i]))
    else:
        for i in range(20):
            colors.append(mpl.colors.rgb2hex(plt.get_cmap("tab20").colors[i]))
        for i in range(20):
            colors.append(mpl.colors.rgb2hex(plt.get_cmap("tab20b").colors[i]))
        for i in range(20):
            colors.append(mpl.colors.rgb2hex(plt.get_cmap("tab20c").colors[i]))
        if k / 60 > 1:
            colors = colors * int(np.ceil(k / 60))

    return colors
