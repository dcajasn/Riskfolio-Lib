""""""  #
"""
Copyright (c) 2020-2024, Dany Cajas
All rights reserved.
This work is licensed under BSD 3-Clause "New" or "Revised" License.
License available at https://github.com/dcajasn/Riskfolio-Lib/blob/master/LICENSE.txt
"""

import numpy as np
import pandas as pd
from riskfolio.external.functions import *
from itertools import product


def duplication_matrix(n: int):
    r"""
    Calculate duplication matrix of size "n" as shown in :cite:`d-Magnus1980`.

    Parameters
    ----------
    n : int
        Number of assets.

    Returns
    -------
    D: np.ndarray
        Duplication matrix
    """
    return cpp_duplication_matrix(n)


def duplication_elimination_matrix(n: int):
    r"""
    Calculate duplication elimination matrix of size "n" as shown in :cite:`d-Magnus1980`.

    Parameters
    ----------
    n : int
        Number of assets.

    Returns
    -------
    L: np.ndarray
        Duplication matrix
    """
    return cpp_duplication_elimination_matrix(n)


def duplication_summation_matrix(n: int):
    r"""
    Calculate duplication summation matrix of size "n" as shown in :cite:`d-Cajas4`.

    Parameters
    ----------
    n : int
        Number of assets.

    Returns
    -------
    S: np.ndarray
        Duplication summation matrix.
    """
    return cpp_duplication_summation_matrix(n)


def commutation_matrix(T: int, n: int):
    r"""
    Calculate commutation matrix of size T x n.

    Parameters
    ----------
    T : int
        Number of rows.

    n : int
        Number of columns.

    Returns
    -------
    K: np.ndarray
        Duplication summation matrix.
    """
    return cpp_commutation_matrix(T, n)


def coskewness_matrix(Y: np.ndarray):
    r"""
    Calculates coskewness rectangular matrix as shown in :cite:`d-Cajas4`.

    Parameters
    ----------
    Y : ndarray or dataframe
        Returns series of shape n_sample x n_features.

    Returns
    -------
    M3 : ndarray
        The lower semi coskewness rectangular matrix.

    Raises
    ------
        ValueError when the value cannot be calculated.

    """
    flag = False
    if isinstance(Y, pd.DataFrame):
        assets = Y.columns.tolist()
        cols = list(product(assets, assets))
        cols = [str(y) + " - " + str(x) for x, y in cols]
        flag = True

    Y_ = np.array(Y, ndmin=2)
    M3 = cpp_coskewness_matrix(Y_, semi=False)

    if flag:
        M3 = pd.DataFrame(M3, index=assets, columns=cols)

    return M3


def semi_coskewness_matrix(Y: np.ndarray):
    r"""
    Calculates lower semi coskewness rectangular matrix as shown in :cite:`d-Cajas4`.

    Parameters
    ----------
    Y : ndarray or dataframe
        Returns series of shape n_samples x n_features.

    Returns
    -------
    s_M3 : ndarray
        The lower semi coskewness rectangular matrix.

    Raises
    ------
        ValueError when the value cannot be calculated.

    """
    flag = False
    if isinstance(Y, pd.DataFrame):
        assets = Y.columns.tolist()
        cols = list(product(assets, assets))
        cols = [str(y) + " - " + str(x) for x, y in cols]
        flag = True

    Y_ = np.array(Y, ndmin=2)
    s_M3 = cpp_coskewness_matrix(Y_, semi=True)

    if flag:
        s_M3 = pd.DataFrame(s_M3, index=assets, columns=cols)

    return s_M3


def cokurtosis_matrix(Y: np.ndarray):
    r"""
    Calculates cokurtosis square matrix as shown in :cite:`d-Cajas4`.

    Parameters
    ----------
    Y : ndarray or dataframe
        Returns series of shape n_samples x n_features.

    Returns
    -------
    S4 : ndarray
        The cokurtosis square matrix.

    Raises
    ------
        ValueError when the value cannot be calculated.

    """
    flag = False
    if isinstance(Y, pd.DataFrame):
        assets = Y.columns.tolist()
        cols = list(product(assets, assets))
        cols = [str(y) + " - " + str(x) for x, y in cols]
        flag = True

    Y_ = np.array(Y, ndmin=2)
    S4 = cpp_cokurtosis_matrix(Y_, semi=False)

    if flag:
        S4 = pd.DataFrame(S4, index=cols, columns=cols)

    return S4


def semi_cokurtosis_matrix(Y):
    r"""
    Calculates lower semi cokurtosis square matrix as shown in :cite:`d-Cajas4`.

    Parameters
    ----------
    Y : ndarray or dataframe
        Returns series of shape n_sample x n_features.

    Returns
    -------
    s_S4 : ndarray
        The lower semi cokurtosis square matrix.

    Raises
    ------
        ValueError when the value cannot be calculated.

    """
    flag = False
    if isinstance(Y, pd.DataFrame):
        assets = Y.columns.tolist()
        cols = list(product(assets, assets))
        cols = [str(y) + " - " + str(x) for x, y in cols]
        flag = True

    Y_ = np.array(Y, ndmin=2)
    s_S4 = cpp_cokurtosis_matrix(Y_, semi=True)

    if flag:
        s_S4 = pd.DataFrame(s_S4, index=cols, columns=cols)

    return s_S4


def k_eigh(Y, k):
    r"""
    Calculates lower semi cokurtosis square matrix as shown in :cite:`d-Cajas4`.

    Parameters
    ----------
    Y : ndarray or dataframe
        Returns series of shape n_sample x n_features.

    Returns
    -------
    s_S4 : ndarray
        The lower semi cokurtosis square matrix.

    Raises
    ------
        ValueError when the value cannot be calculated.
    """
    Y_ = np.array(Y, ndmin=2)
    eigvalues, eigvectors = cpp_k_eigh(Y_, k)

    return eigvalues, eigvectors


def d_corr(X, Y):
    r"""
    Calculates the distance correlation of X and Y.

    Parameters
    ----------
    X : ndarray or dataframe
        Returns series of shape n_sample x n_features.
    Y : ndarray or dataframe
        Returns series of shape n_sample x n_features.

    Returns
    -------
    value : float
        Distance correlation.

    Raises
    ------
        ValueError when the value cannot be calculated.
    """
    X_ = np.array(X, ndmin=2)
    Y_ = np.array(Y, ndmin=2)
    value = cpp_dcorr(X_, Y_)

    return value


def d_corr_matrix(Y):
    r"""
    Calculates the distance correlation matrix of matrix of variables Y.

    Parameters
    ----------
    Y : ndarray or dataframe
        Returns series of shape n_sample x n_features.

    Returns
    -------
    value : float
        Distance correlation.

    Raises
    ------
        ValueError when the value cannot be calculated.

    """
    Y_ = np.array(Y, ndmin=2)
    value = cpp_dcorr_matrix(Y_)

    return value
