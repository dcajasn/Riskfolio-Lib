""""""  #
"""
Copyright (c) 2020-2023, Dany Cajas
All rights reserved.
This work is licensed under BSD 3-Clause "New" or "Revised" License.
License available at https://github.com/dcajasn/Riskfolio-Lib/blob/master/LICENSE.txt
"""

from riskfolio.external.functions import *


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
