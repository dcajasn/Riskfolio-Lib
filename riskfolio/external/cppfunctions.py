from .functions import *


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
    """
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
