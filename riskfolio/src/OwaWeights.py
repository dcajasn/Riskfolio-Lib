""""""  #
"""
Copyright (c) 2020-2024, Dany Cajas
All rights reserved.
This work is licensed under BSD 3-Clause "New" or "Revised" License.
License available at https://github.com/dcajasn/Riskfolio-Lib/blob/master/LICENSE.txt
"""

import numpy as np
import cvxpy as cp
import math
from scipy.special import binom

__all__ = [
    "owa_l_moment",
    "owa_gmd",
    "owa_cvar",
    "owa_wcvar",
    "owa_tg",
    "owa_wr",
    "owa_rg",
    "owa_cvrg",
    "owa_wcvrg",
    "owa_tgrg",
    "owa_l_moment_crm",
]


def owa_l_moment(T, k=2):
    r"""
    Calculate the OWA weights to calculate the kth linear moment (l-moment)
    of a returns series as shown in :cite:`d-Cajas6`.

    Parameters
    ----------
    T : int
        Number of observations of the returns series.
    k : int
        Order of the l-moment. Must be an integer higher or equal than 1.

    Returns
    -------
    value : 1d-array
        An OWA weights vector of size Tx1.
    """
    w = []
    for i in range(1, T + 1):
        a = 0
        for j in range(k):
            a += (-1) ** j * binom(k - 1, j) * binom(i - 1, k - 1 - j) * binom(T - i, j)
        a *= 1 / (k * binom(T, k))
        w.append(a)
    return np.array(w).reshape(-1, 1)


def owa_gmd(T):
    r"""
    Calculate the OWA weights to calculate the Gini mean difference (GMD)
    of a returns series as shown in :cite:`d-Cajas3`.

    Parameters
    ----------
    T : int
        Number of observations of the returns series.

    Returns
    -------
    value : 1d-array
        An OWA weights vector of size Tx1.
    """

    w_ = []
    for i in range(1, T + 1):
        w_.append(2 * i - 1 - T)
    w_ = 2 * np.array(w_) / (T * (T - 1))
    w_ = w_.reshape(-1, 1)

    return w_


def owa_cvar(T, alpha=0.05):
    r"""
    Calculate the OWA weights to calculate the Conditional Value at Risk (CVaR)
    of a returns series as shown in :cite:`d-Cajas3`.

    Parameters
    ----------
    T : int
        Number of observations of the returns series.
    alpha : float, optional
        Significance level of CVaR. The default is 0.05.

    Returns
    -------
    value : 1d-array
        An OWA weights vector of size Tx1.
    """

    k = int(np.ceil(T * alpha)) - 1
    w_ = np.zeros((T, 1))
    w_[:k, :] = -1 / (T * alpha)
    w_[k, :] = -1 - np.sum(w_[:k, :])

    return w_


def owa_wcvar(T, alphas, weights):
    r"""
    Calculate the OWA weights to calculate the Weighted Conditional Value at
    Risk (WCVaR) of a returns series as shown in :cite:`d-Cajas3`.

    Parameters
    ----------
    T : int
        Number of observations of the returns series.
    alphas : list
        List of significance levels of each CVaR model.
    weights : list
        List of weights of each CVaR model.

    Returns
    -------
    value : 1d-array
        An OWA weights vector of size Tx1.
    """

    w_ = 0
    for i, j in zip(alphas, weights):
        w_ += owa_cvar(T, i) * j

    return w_


def owa_tg(T, alpha=0.05, a_sim=100):
    r"""
    Calculate the OWA weights to calculate the Tail Gini of a
    returns series as shown in :cite:`d-Cajas3`.

    Parameters
    ----------
    T : int
        Number of observations of the returns series.
    alpha : float, optional
        Significance level of TaiL Gini. The default is 0.05.
    a_sim : float, optional
        Number of CVaRs used to approximate the Tail Gini. The default is 100.

    Returns
    -------
    value : 1d-array
        A OWA weights vector of size Tx1.
    """

    alphas = np.linspace(alpha, 0.0001, a_sim)[::-1]
    w_ = [(alphas[1] - 0) * alphas[0] / alphas[-1] ** 2]
    for i in range(1, len(alphas) - 1):
        w_.append((alphas[i + 1] - alphas[i - 1]) * alphas[i] / alphas[-1] ** 2)
    w_.append((alphas[-1] - alphas[-2]) / alphas[-1])
    w_ = owa_wcvar(T, alphas, w_)

    return w_


def owa_wr(T):
    r"""
    Calculate the OWA weights to calculate the Worst realization (minimum)
    of a returns series as shown in :cite:`d-Cajas3`.

    Parameters
    ----------
    T : int
        Number of observations of the returns series.

    Returns
    -------
    value : 1d-array
        A OWA weights vector of size Tx1.
    """

    w_ = np.zeros((T, 1))
    w_[0, :] = -1

    return w_


def owa_rg(T):
    r"""
    Calculate the OWA weights to calculate the range of a returns series
    as shown in :cite:`d-Cajas3`.

    Parameters
    ----------
    T : int
        Number of observations of the returns series.

    Returns
    -------
    value : 1d-array
        A OWA weights vector of size Tx1.
    """

    w_ = np.zeros((T, 1))
    w_[0, :] = -1
    w_[-1, :] = 1

    return w_


def owa_cvrg(T, alpha=0.05, beta=None):
    r"""
    Calculate the OWA weights to calculate the CVaR range of a returns series
    as shown in :cite:`d-Cajas3`.

    Parameters
    ----------
    T : int
        Number of observations of the returns series.
    alpha : float, optional
        Significance level of CVaR of losses. The default is 0.05.
    beta : float, optional
        Significance level of CVaR of gains. If None it duplicates alpha.
        The default is None.

    Returns
    -------
    value : 1d-array
        A OWA weights vector of size Tx1.
    """

    if beta is None:
        beta = alpha

    w_ = owa_cvar(T, alpha) - owa_cvar(T, beta)[::-1]

    return w_


def owa_wcvrg(T, alphas, weights_a, betas=None, weights_b=None):
    r"""
    Calculate the OWA weights to calculate the WCVaR range of a returns series
    as shown in :cite:`d-Cajas3`.

    Parameters
    ----------
    T : int
        Number of observations of the returns series.
    alphas : list
        List of significance levels of each CVaR of losses model.
    weights_a : list
        List of weights of each CVaR of losses model.
    betas : list, optional
        List of significance levels of each CVaR of gains model. If None it duplicates alpha.
        The default is None.
    weights_b : list, optional
        List of weights of each CVaR of gains model. If None it duplicates weights_a.
        The default is None.

    Returns
    -------
    value : 1d-array
        A OWA weights vector of size Tx1.
    """

    if betas is None or weights_b is None:
        betas = alphas
        weights_b = weights_a

    w_ = owa_wcvar(T, alphas, weights_a) - owa_wcvar(T, betas, weights_b)[::-1]

    return w_


def owa_tgrg(T, alpha=0.05, a_sim=100, beta=None, b_sim=None):
    r"""
    Calculate the OWA weights to calculate the Tail Gini range of a returns
    series as shown in :cite:`d-Cajas3`.

    Parameters
    ----------
    T : int
        Number of observations of the returns series.
    alpha : float, optional
        Significance level of Tail Gini of losses. The default is 0.05.
    a_sim : float, optional
        Number of CVaRs used to approximate Tail Gini of losses. The default is 100.
    beta : float, optional
        Significance level of Tail Gini of gains. If None it duplicates alpha value.
        The default is None.
    b_sim : float, optional
        Number of CVaRs used to approximate Tail Gini of gains. If None it duplicates a_sim value.
        The default is None.

    Returns
    -------
    value : 1d-array
        A OWA weights vector of size Tx1.
    """

    if beta is None:
        beta = alpha
    if b_sim is None:
        b_sim = a_sim

    w_ = owa_tg(T, alpha, a_sim) - owa_tg(T, beta, b_sim)[::-1]

    return w_


def owa_l_moment_crm(T, k=4, method="MSD", g=0.5, max_phi=0.5, solver='CLARABEL'):
    r"""
    Calculate the OWA weights to calculate a convex risk measure that considers
    higher linear moments or L-moments as shown in :cite:`d-Cajas6`.

    Parameters
    ----------
    T : int
        Number of observations of the returns series.
    k : int
        Order of the l-moment. Must be an integer higher or equal than 2.
    method : str, optional
        Method to calculate the weights used to combine the l-moments with order higher than 2.
        The default value is 'MSD'. Possible values are:

        - 'CRRA': Normalized Constant Relative Risk Aversion coefficients.
        - 'ME': Maximum Entropy.
        - 'MSS': Minimum Sum Squares.
        - 'MSD': Minimum Square Distance.

    g : float, optional
        Risk aversion coefficient of CRRA utility function. The default is 0.5.
    max_phi : float, optional
        Maximum weight constraint of L-moments.
        The default is 0.5.
    solver: str, optional
        Solver available for CVXPY. Used to calculate 'ME', 'MSS' and 'MSD' weights.
        The default value is 'CLARABEL'.

    Returns
    -------
    value : 1d-array
        A OWA weights vector of size Tx1.
    """

    if k < 2 or (not isinstance(k, int)):
        raise ValueError("k must be an integer higher equal than 2")
    if method not in ["CRRA", "ME", "MSS", "MSD"]:
        raise ValueError("Available methods are 'CRRA', 'ME', 'MSS' and 'MSD'")
    if g >= 1 or g <= 0:
        raise ValueError("The risk aversion coefficient mus be between 0 and 1")
    if max_phi >= 1 or max_phi <= 0:
        raise ValueError(
            "The constraint on maximum weight of L-moments must be between 0 and 1"
        )

    ws = np.empty((T, 0))
    for i in range(2, k + 1):
        w_i = (-1) ** i * owa_l_moment(T, i)
        ws = np.concatenate([ws, w_i], axis=1)

    if method == "CRRA":
        phis = []
        e = 1
        for i in range(1, k):
            e *= g + i - 1
            phis.append(e / math.factorial(i + 1))
        phis = np.array(phis)
        phis = phis / np.sum(phis)
        phis = phis.reshape(-1, 1)
        a = ws @ phis

        w = np.zeros_like(a)
        w[0] = a[0]
        for i in range(1, len(a)):
            w[i, 0] = np.max(a[: i + 1, 0])

    else:
        theta = cp.Variable((T, 1))
        n = ws.shape[1]
        phi = cp.Variable((n, 1))

        constraints = [
            cp.sum(phi) == 1,
            theta == ws @ phi,
            phi <= max_phi,
            phi >= 0,
            phi[1:] <= phi[:-1],
            theta[1:] >= theta[:-1],
        ]

        if method == "ME":
            theta_ = cp.Variable((T, 1))
            obj = cp.sum(cp.entr(theta_)) * 1000
            constraints += [
                theta_ >= theta,
                theta_ >= -theta,
            ]
            objective = cp.Maximize(obj)
        elif method == "MSS":
            obj = cp.pnorm(theta, p=2) * 1000
            objective = cp.Minimize(obj)
        elif method == "MSD":
            obj = cp.pnorm(theta[1:] - theta[:-1], p=2) * 1000
            objective = cp.Minimize(obj)

        problem = cp.Problem(objective, constraints)
        if solver is not None:
            problem.solve(solver=solver)
        else:
            problem.solve()

        phis = phi.value
        phis = phis / np.sum(phis)
        phis = phis.reshape(-1, 1)
        w = ws @ phis

    return w
