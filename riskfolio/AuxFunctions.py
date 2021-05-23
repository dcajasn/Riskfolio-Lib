import numpy as np
import pandas as pd
from scipy import linalg as LA
from statsmodels.stats.correlation_tools import cov_nearest
from scipy.sparse import csr_matrix
from scipy.spatial.distance import pdist, squareform
import scipy.cluster.hierarchy as hr

###############################################################################
# Some Aditional Functions
###############################################################################


def is_pos_def(cov, threshold=1e-8):
    r"""
    Indicate if a matrix is positive (semi)definite.

    Parameters
    ----------
    cov : nd-array of shape (n_features, n_features)
        Features covariance matrix, where n_features is the number of features.

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
    cov : nd-array of shape (n_features, n_features)
        Assets covariance matrix, where n_features is the number of features.

    Returns
    -------
    corr : nd-array
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
    cov : nd-array of shape (n_features, n_features)
        Features covariance matrix, where n_features is the number of features.
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
    cov : nd-array of shape (n_features, n_features)
        Assets covariance matrix, where n_features is the number of features.

    Returns
    -------
    a : nd-array
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
    cov : nd-array of shape (n_features, n_features)
        Assets covariance matrix, where n_features is the number of features.

    Returns
    -------
    K : nd-array
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


def dcorr(X, Y):
    r"""
    Calculate the distance correlation between two variables :cite:`d-Szekely`.

    Parameters
    ----------
    X : 1d-array of shape (n_sample, 1)
        Returns series, must have Tx1 size.
    Y : 1d-array of shape (n_sample, 1)
        Returns series, must have Tx1 size.

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
    X : 2d-array of shape (n_sample, n_features)
        Returns series, must have Txn size.
    Returns
    -------
    corr : 2d-array of shape (n_features, n_features)
        The distance correlation matrix of n variables.

    Raises
    ------
        ValueError when the value cannot be calculated.

    """
    flag = False
    if isinstance(X, pd.DataFrame):
        cols = X.columns.tolist()
        flag = True

    n = X.shape[1]
    corr = np.ones((n, n))
    indices = np.triu_indices(n, 1)

    for i, j in zip(indices[0], indices[1]):
        corr[i, j] = dcorr(X.iloc[:, i], X.iloc[:, j])
        corr[j, i] = corr[i, j]

    if flag:
        corr = pd.DataFrame(corr, index=cols, columns=cols)
    else:
        corr = pd.DataFrame(corr)

    return corr


def two_diff_gap_stat(corr, dist, clusters, max_k=10):
    r"""
    Calculate the optimal number of clusters based on the two difference gap
    statistic :cite:`d-twogap`.

    Parameters
    ----------
    corr : DataFrame
        A correlation matrix.
    dist : str, optional
        A distance measure based on the correlation matrix.
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
    cluster_lvls = pd.DataFrame(hr.cut_tree(clusters), index=corr.columns)
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
            D = np.nan_to_num(cluster_pdist.mean())
            D_list.append(D)  # append to list

        W_k = np.sum(D_list)
        W_list.append(W_k)

    W_list = pd.Series(W_list)
    n = corr.shape[0]
    limit_k = int(min(max_k, np.sqrt(n)))
    gaps = W_list.shift(2) + W_list - 2 * W_list.shift(1)
    gaps = gaps[0:limit_k]
    k = int(gaps.idxmax() + 2)

    return k
