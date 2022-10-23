""""""  #
"""
This code is mainly based on Yinsen Miao's work available in:
https://github.com/yinsenm/gerber/blob/af04c2ee5adf342393b028b85ab5546f31c0c8d3/src/gerber.py
"""

import numpy as np
import pandas as pd
import riskfolio.src.AuxFunctions as af


def gerber_cov_stat0(X, threshold=0.5):
    r"""
    Compute Gerber covariance Statistics 0 or original Gerber statistics
    :cite: `d-Gerber2021`, not always PSD, however this function fixes the
    covariance matrix finding the nearest covariance matrix that is positive
    semidefinite.

    Parameters
    ----------
    X : ndarray
        Returns series of shape n_sample x n_features.

    threshold : float
        threshold: threshold is between 0 and 1

    Returns
    -------
    value : bool
        Gerber covariance matrix of shape (n_features, n_features), where
        n_features is the number of features.

    Raises
    ------
        ValueError when the value cannot be calculated.

    """
    # Threshold shall between 0 and 1
    assert 1 > threshold > 0

    flag = False
    if isinstance(X, pd.DataFrame):
        cols = X.columns.tolist()
        X1 = X.to_numpy()
        flag = True
    else:
        X1 = X.copy()

    n, p = X1.shape
    sd_vec = X1.std(axis=0).reshape((p, 1))
    cov = np.zeros((p, p))  # Store correlation matrix
    corr = np.zeros((p, p))  # Store correlation matrix

    for i in range(p):
        for j in range(i + 1):
            neg = 0
            pos = 0
            for k in range(n):
                if (
                    (X1[k, i] >= threshold * sd_vec[i])
                    and (X1[k, j] >= threshold * sd_vec[j])
                ) or (
                    (X1[k, i] <= -threshold * sd_vec[i])
                    and (X1[k, j] <= -threshold * sd_vec[j])
                ):
                    pos += 1
                elif (
                    (X1[k, i] >= threshold * sd_vec[i])
                    and (X1[k, j] <= -threshold * sd_vec[j])
                ) or (
                    (X1[k, i] <= -threshold * sd_vec[i])
                    and (X1[k, j] >= threshold * sd_vec[j])
                ):
                    neg += 1

            # Compute Gerber correlation matrix
            corr[i, j] = (pos - neg) / (pos + neg)
            corr[j, i] = corr[i, j]

    cov = corr * np.outer(sd_vec, sd_vec)
    if af.is_pos_def(cov) == False:
        cov = af.cov_fix(cov, method="clipped")

    if flag:
        cov = pd.DataFrame(cov, index=cols, columns=cols)

    return cov


def gerber_cov_stat1(X, threshold=0.5):
    r"""
    Compute Gerber covariance Statistics 1 :cite: `d-Gerber2021`.

    Parameters
    ----------
    X : ndarray
        Returns series of shape n_sample x n_features.

    threshold : float
        threshold: threshold is between 0 and 1

    Returns
    -------
    value : bool
        Gerber covariance matrix of shape (n_features, n_features), where
        n_features is the number of features.

    Raises
    ------
        ValueError when the value cannot be calculated.

    """
    # Threshold shall between 0 and 1
    assert 1 > threshold > 0

    flag = False
    if isinstance(X, pd.DataFrame):
        cols = X.columns.tolist()
        X1 = X.to_numpy()
        flag = True
    else:
        X1 = X.copy()

    n, p = X1.shape
    sd_vec = X1.std(axis=0).reshape((p, 1))
    corr = np.zeros((p, p))  # Store correlation matrix

    for i in range(p):
        for j in range(i + 1):
            neg = 0
            pos = 0
            nn = 0
            for k in range(n):
                if (
                    (X1[k, i] >= threshold * sd_vec[i])
                    and (X1[k, j] >= threshold * sd_vec[j])
                ) or (
                    (X1[k, i] <= -threshold * sd_vec[i])
                    and (X1[k, j] <= -threshold * sd_vec[j])
                ):
                    pos += 1
                elif (
                    (X1[k, i] >= threshold * sd_vec[i])
                    and (X1[k, j] <= -threshold * sd_vec[j])
                ) or (
                    (X1[k, i] <= -threshold * sd_vec[i])
                    and (X1[k, j] >= threshold * sd_vec[j])
                ):
                    neg += 1
                elif (
                    abs(X1[k, i]) < threshold * sd_vec[i]
                    and abs(X1[k, j]) < threshold * sd_vec[j]
                ):
                    nn += 1

            # Compute Gerber correlation matrix
            corr[i, j] = (pos - neg) / (n - nn)
            corr[j, i] = corr[i, j]

    cov = corr * np.outer(sd_vec, sd_vec)

    if flag:
        cov = pd.DataFrame(cov, index=cols, columns=cols)

    return cov


def gerber_cov_stat2(X, threshold=0.5):
    r"""
    Compute Gerber covariance Statistics 2 :cite: `d-Gerber2021`.

    Parameters
    ----------
    X : ndarray
        Returns series of shape n_sample x n_features.

    threshold : float
        threshold: threshold is between 0 and 1

    Returns
    -------
    value : bool
        Gerber covariance mtrix of shape (n_features, n_features), where
        n_features is the number of features.

    Raises
    ------
        ValueError when the value cannot be calculated.

    """
    # Threshold shall between 0 and 1
    assert 1 > threshold > 0

    flag = False
    if isinstance(X, pd.DataFrame):
        cols = X.columns.tolist()
        X1 = X.to_numpy()
        flag = True
    else:
        X1 = X.copy()

    n, p = X1.shape
    sd_vec = X1.std(axis=0).reshape((p, 1))
    U = np.copy(X1)
    D = np.copy(X1)

    # Update U and D matrix
    for i in range(p):
        U[:, i] = U[:, i] >= sd_vec[i] * threshold
        D[:, i] = D[:, i] <= -sd_vec[i] * threshold

    # Update concordant matrix
    N_CONC = U.transpose() @ U + D.transpose() @ D

    # Update discordant matrix
    N_DISC = U.transpose() @ D + D.transpose() @ U
    H = N_CONC - N_DISC
    h = np.sqrt(H.diagonal())
    h = h.reshape((p, 1))

    corr = H / (h @ h.transpose())

    cov = corr * np.outer(sd_vec, sd_vec)

    if flag:
        cov = pd.DataFrame(cov, index=cols, columns=cols)

    return cov
