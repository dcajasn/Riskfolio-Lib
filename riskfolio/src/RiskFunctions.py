""""""  #

"""
Copyright (c) 2020-2025, Dany Cajas
All rights reserved.
This work is licensed under BSD 3-Clause "New" or "Revised" License.
License available at https://github.com/dcajasn/Riskfolio-Lib/blob/master/LICENSE.txt
"""

import numpy as np
import pandas as pd
import cvxpy as cp
import riskfolio.src.OwaWeights as owa
import riskfolio.src.ParamsEstimation as pe
from scipy.optimize import minimize
from scipy.optimize import Bounds
from scipy.linalg import null_space
from numpy.linalg import pinv
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import warnings


__all__ = [
    "MAD",
    "SemiDeviation",
    "Kurtosis",
    "SemiKurtosis",
    "VaR_Hist",
    "CVaR_Hist",
    "WR",
    "LPM",
    "Entropic_RM",
    "EVaR_Hist",
    "RLVaR_Hist",
    "MDD_Abs",
    "ADD_Abs",
    "DaR_Abs",
    "CDaR_Abs",
    "EDaR_Abs",
    "RLDaR_Abs",
    "UCI_Abs",
    "MDD_Rel",
    "ADD_Rel",
    "DaR_Rel",
    "CDaR_Rel",
    "EDaR_Rel",
    "RLDaR_Rel",
    "UCI_Rel",
    "GMD",
    "TG",
    "RG",
    "VRG",
    "CVRG",
    "TGRG",
    "EVRG",
    "RVRG",
    "L_Moment",
    "L_Moment_CRM",
    "NEA",
    "Sharpe_Risk",
    "Sharpe",
    "Risk_Contribution",
    "Risk_Margin",
    "Factors_Risk_Contribution",
    "BrinsonAttribution",
]


def MAD(X):
    r"""
    Calculate the Mean Absolute Deviation (MAD) of a returns series.

    .. math::
        \text{MAD}(X) = \frac{1}{T}\sum_{t=1}^{T}
        | X_{t} - \mathbb{E}(X_{t}) |

    Parameters
    ----------
    X : 1d-array
        Returns series, must have Tx1 size.

    Returns
    -------
    value : float
        MAD of a returns series.

    """

    a = np.array(X, ndmin=2)
    if a.shape[0] == 1 and a.shape[1] > 1:
        a = a.T
    if a.shape[0] > 1 and a.shape[1] > 1:
        raise ValueError("returns must have Tx1 size")

    T, N = a.shape
    mu = np.mean(a, axis=0).reshape(1, -1)
    mu = np.repeat(mu, T, axis=0)
    value = a - mu
    value = np.mean(np.absolute(value), axis=0)
    value = np.array(value).item()

    return value


def SemiDeviation(X):
    r"""
    Calculate the Semi Deviation of a returns series.

    .. math::
        \text{SemiDev}(X) = \left [ \frac{1}{T-1}\sum_{t=1}^{T}
        \min (X_{t} - \mathbb{E}(X_{t}), 0)^2 \right ]^{1/2}

    Parameters
    ----------
    X : 1d-array
        Returns series, must have Tx1 size.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    Returns
    -------
    value : float
        Semi Deviation of a returns series.
    """

    a = np.array(X, ndmin=2)
    if a.shape[0] == 1 and a.shape[1] > 1:
        a = a.T
    if a.shape[0] > 1 and a.shape[1] > 1:
        raise ValueError("returns must have Tx1 size")

    T, N = a.shape
    mu = np.mean(a, axis=0).reshape(1, -1)
    mu = np.repeat(mu, T, axis=0)
    value = mu - a
    value = np.sum(np.power(value[np.where(value >= 0)], 2)) / (T - 1)
    value = np.power(value, 0.5).item()

    return value


def Kurtosis(X):
    r"""
    Calculate the Square Root Kurtosis of a returns series.

    .. math::
        \text{Kurt}(X) = \left [ \frac{1}{T}\sum_{t=1}^{T}
        (X_{t} - \mathbb{E}(X_{t}))^{4} \right ]^{1/2}

    Parameters
    ----------
    X : 1d-array
        Returns series, must have Tx1 size.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    Returns
    -------
    value : float
        Square Root Kurtosis of a returns series.
    """

    a = np.array(X, ndmin=2)
    if a.shape[0] == 1 and a.shape[1] > 1:
        a = a.T
    if a.shape[0] > 1 and a.shape[1] > 1:
        raise ValueError("returns must have Tx1 size")

    T, N = a.shape
    mu = np.mean(a, axis=0).reshape(1, -1)
    mu = np.repeat(mu, T, axis=0)
    value = mu - a
    value = np.sum(np.power(value, 4)) / T
    value = np.power(value, 0.5).item()

    return value


def SemiKurtosis(X):
    r"""
    Calculate the Semi Square Root Kurtosis of a returns series.

    .. math::
        \text{SemiKurt}(X) = \left [ \frac{1}{T}\sum_{t=1}^{T}
        \min (X_{t} - \mathbb{E}(X_{t}), 0)^{4} \right ]^{1/2}

    Parameters
    ----------
    X : 1d-array
        Returns series, must have Tx1 size.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    Returns
    -------
    value : float
        Semi Square Root Kurtosis of a returns series.
    """

    a = np.array(X, ndmin=2)
    if a.shape[0] == 1 and a.shape[1] > 1:
        a = a.T
    if a.shape[0] > 1 and a.shape[1] > 1:
        raise ValueError("returns must have Tx1 size")

    T, N = a.shape
    mu = np.mean(a, axis=0).reshape(1, -1)
    mu = np.repeat(mu, T, axis=0)
    value = mu - a
    value = np.sum(np.power(value[np.where(value >= 0)], 4)) / T
    value = np.power(value, 0.5).item()

    return value


def VaR_Hist(X, alpha=0.05):
    r"""
    Calculate the Value at Risk (VaR) of a returns series.

    .. math::
        \text{VaR}_{\alpha}(X) = -\inf_{t \in (0,T)} \left \{ X_{t} \in
        \mathbb{R}: F_{X}(X_{t})>\alpha \right \}

    Parameters
    ----------
    X : 1d-array
        Returns series, must have Tx1 size.
    alpha : float, optional
        Significance level of VaR. The default is 0.05.
    Raises
    ------
    ValueError
        When the value cannot be calculated.

    Returns
    -------
    value : float
        VaR of a returns series.
    """

    a = np.array(X, ndmin=2)
    if a.shape[0] == 1 and a.shape[1] > 1:
        a = a.T
    if a.shape[0] > 1 and a.shape[1] > 1:
        raise ValueError("returns must have Tx1 size")

    sorted_a = np.sort(a, axis=0)
    index = int(np.ceil(alpha * len(sorted_a)) - 1)
    value = -sorted_a[index]
    value = np.array(value).item()

    return value


def CVaR_Hist(X, alpha=0.05):
    r"""
    Calculate the Conditional Value at Risk (CVaR) of a returns series.

    .. math::
        \text{CVaR}_{\alpha}(X) = \text{VaR}_{\alpha}(X) +
        \frac{1}{\alpha T} \sum_{t=1}^{T} \max(-X_{t} -
        \text{VaR}_{\alpha}(X), 0)

    Parameters
    ----------
    X : 1d-array
        Returns series, must have Tx1 size.
    alpha : float, optional
        Significance level of CVaR. The default is 0.05.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    Returns
    -------
    value : float
        CVaR of a returns series.
    """

    a = np.array(X, ndmin=2)
    if a.shape[0] == 1 and a.shape[1] > 1:
        a = a.T
    if a.shape[0] > 1 and a.shape[1] > 1:
        raise ValueError("returns must have Tx1 size")

    sorted_a = np.sort(a, axis=0)
    index = int(np.ceil(alpha * len(sorted_a)) - 1)
    sum_var = 0
    for i in range(0, index + 1):
        sum_var = sum_var + sorted_a[i] - sorted_a[index]

    value = -sorted_a[index] - sum_var / (alpha * len(sorted_a))
    value = np.array(value).item()

    return value


def WR(X):
    r"""
    Calculate the Worst Realization (WR) or Worst Scenario of a returns series.

    .. math::
        \text{WR}(X) = \max(-X)

    Parameters
    ----------
    X : 1d-array
        Returns series, must have Tx1 size.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    Returns
    -------
    value : float
        WR of a returns series.

    """

    a = np.array(X, ndmin=2)
    if a.shape[0] == 1 and a.shape[1] > 1:
        a = a.T
    if a.shape[0] > 1 and a.shape[1] > 1:
        raise ValueError("returns must have Tx1 size")

    sorted_a = np.sort(a, axis=0)
    value = -sorted_a[0]
    value = np.array(value).item()

    return value


def LPM(X, MAR=0, p=1):
    r"""
    Calculate the First or Second Lower Partial Moment of a returns series.

    .. math::
        \text{LPM}(X, \text{MAR}, 1) &= \frac{1}{T}\sum_{t=1}^{T}
        \max(\text{MAR} - X_{t}, 0) \\
        \text{LPM}(X, \text{MAR}, 2) &= \left [ \frac{1}{T-1}\sum_{t=1}^{T}
        \max(\text{MAR} - X_{t}, 0)^{2} \right ]^{\frac{1}{2}} \\


    Where:

    :math:`\text{MAR}` is the minimum acceptable return.
    :math:`p` is the order of the :math:`\text{LPM}`.

    Parameters
    ----------
    X : 1d-array
        Returns series, must have Tx1 size.
    MAR : float, optional
        Minimum acceptable return. The default is 0.
    p : float, optional can be {1,2}
        order of the :math:`\text{LPM}`. The default is 1.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    Returns
    -------
    value : float
        p-th Lower Partial Moment of a returns series.

    """

    a = np.array(X, ndmin=2)
    if a.shape[0] == 1 and a.shape[1] > 1:
        a = a.T
    if a.shape[0] > 1 and a.shape[1] > 1:
        raise ValueError("returns must have Tx1 size")
    if p not in [1, 2]:
        raise ValueError("p can only be 1 or 2")

    value = MAR - a

    if p == 2:
        n = value.shape[0] - 1
    else:
        n = value.shape[0]

    value = np.sum(np.power(value[np.where(value >= 0)], p)) / n
    value = np.power(value, 1 / p).item()

    return value


def Entropic_RM(X, z=1, alpha=0.05):
    r"""
    Calculate the Entropic Risk Measure (ERM) of a returns series.

    .. math::
        \text{ERM}_{\alpha}(X) = z\ln \left (\frac{M_X(z^{-1})}{\alpha} \right )

    Where:

    :math:`M_X(z)` is the moment generating function of X.

    Parameters
    ----------
    X : 1d-array
        Returns series, must have Tx1 size.
    z : float, optional
        Risk aversion parameter, must be greater than zero. The default is 1.
    alpha : float, optional
        Significance level of EVaR. The default is 0.05.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    Returns
    -------
    value : float
        ERM of a returns series.

    """

    a = np.array(X, ndmin=2)
    if a.shape[0] == 1 and a.shape[1] > 1:
        a = a.T
    if a.shape[0] > 1 and a.shape[1] > 1:
        raise ValueError("returns must have Tx1 size")

    value = np.mean(np.exp(-1 / z * a), axis=0)
    value = z * (np.log(value) + np.log(1 / alpha))
    value = np.array(value).item()

    return value


def _Entropic_RM(z, X, alpha=0.05):
    a = np.array(X, ndmin=2)
    if a.shape[0] == 1 and a.shape[1] > 1:
        a = a.T
    if a.shape[0] > 1 and a.shape[1] > 1:
        raise ValueError("returns must have Tx1 size")

    a = a.flatten()
    value = np.mean(np.exp(-1 / z * a), axis=0)
    value = z * (np.log(value) + np.log(1 / alpha))
    value = np.array(value).item()

    return value


def EVaR_Hist(X, alpha=0.05, solver="CLARABEL"):
    r"""
    Calculate the Entropic Value at Risk (EVaR) of a returns series.

    .. math::
        \text{EVaR}_{\alpha}(X) = \inf_{z>0} \left \{ z
        \ln \left (\frac{M_X(z^{-1})}{\alpha} \right ) \right \}

    Where:

    :math:`M_X(t)` is the moment generating function of X.

    Parameters
    ----------
    X : 1d-array
        Returns series, must have Tx1 size.
    alpha : float, optional
        Significance level of EVaR. The default is 0.05.
    solver: str, optional
        Solver available for CVXPY that supports exponential cone programming.
        Used to calculate EVaR, EVRG and EDaR. The default value is 'CLARABEL'.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    Returns
    -------
    (value, z) : tuple
        EVaR of a returns series and value of z that minimize EVaR.

    """

    a = np.array(X, ndmin=2)
    if a.shape[0] == 1 and a.shape[1] > 1:
        a = a.T
    if a.shape[0] > 1 and a.shape[1] > 1:
        raise ValueError("returns must have Tx1 size")

    T, N = a.shape

    # Primal Formulation
    t = cp.Variable((1, 1))
    z = cp.Variable((1, 1), nonneg=True)
    ui = cp.Variable((T, 1))
    ones = np.ones((T, 1))

    constraints = [
        cp.sum(ui) <= z,
        cp.ExpCone(-a * 1000 - t * 1000, ones @ z * 1000, ui * 1000),
    ]

    risk = t + z * np.log(1 / (alpha * T))
    objective = cp.Minimize(risk * 1000)
    prob = cp.Problem(objective, constraints)

    try:
        if solver in ["CLARABEL", "MOSEK", "SCS"]:
            prob.solve(solver=solver)
        else:
            prob.solve()
    except:
        pass

    if risk.value is None:
        value = None
    else:
        value = risk.value.item()
        t = z.value

    if value is None:
        warnings.filterwarnings("ignore")

        bnd = Bounds([-1e-24], [np.inf])
        result = minimize(
            _Entropic_RM, [1], args=(X, alpha), method="SLSQP", bounds=bnd, tol=1e-12
        )
        t = result.x
        t = t.item()
        value = _Entropic_RM(t, X, alpha)

    return (value, t)


def RLVaR_Hist(X, alpha=0.05, kappa=0.3, solver="CLARABEL"):
    r"""
    Calculate the Relativistic Value at Risk (RLVaR) of a returns series.
    I recommend only use this function with MOSEK solver.

    .. math::
        \text{RLVaR}^{\kappa}_{\alpha}(X) & = \left \{
        \begin{array}{ll}
        \underset{z, t, \psi, \theta,  \varepsilon, \omega}{\text{inf}} & t + z \ln_{\kappa} \left ( \frac{1}{\alpha T} \right ) + \sum^T_{i=1} \left ( \psi_{i} + \theta_{i}  \right ) \\
        \text{s.t.} & -X  - t + \varepsilon + \omega \leq 0\\
        & z \geq 0 \\
        & \left ( z \left ( \frac{1+\kappa}{2\kappa} \right ), \psi_{i} \left ( \frac{1+\kappa}{\kappa} \right ), \varepsilon_{i} \right) \in \mathcal{P}_3^{1/(1+\kappa),\, \kappa/(1+\kappa)} \\
        & \left ( \omega_{i}\left ( \frac{1}{1-\kappa} \right ), \theta_{i}\left ( \frac{1}{\kappa} \right),  -z \left ( \frac{1}{2\kappa} \right ) \right ) \in \mathcal{P}_3^{1-\kappa,\, \kappa} \\
        \end{array} \right .

    Where:

    :math:`\mathcal{P}_3^{\alpha,\, 1-\alpha}` is the power cone 3D.

    :math:`\kappa` is the deformation parameter.

    Parameters
    ----------
    X : 1d-array
        Returns series, must have Tx1 size.
    alpha : float, optional
        Significance level of EVaR. The default is 0.05.
    kappa : float, optional
        Deformation parameter of RLVaR, must be between 0 and 1. The default is 0.3.
    solver: str, optional
        Solver available for CVXPY that supports power cone programming. Used
        to calculate RLVaR and RLDaR. The default value is 'CLARABEL'.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    Returns
    -------
    value : tuple
        RLVaR of a returns series.

    """

    a = np.array(X * 100, ndmin=2)
    if a.shape[0] == 1 and a.shape[1] > 1:
        a = a.T
    if a.shape[0] > 1 and a.shape[1] > 1:
        raise ValueError("returns must have Tx1 size")

    T, N = a.shape

    # Dual Formulation
    Z = cp.Variable((T, 1))
    nu = cp.Variable((T, 1))
    tau = cp.Variable((T, 1))
    ones = np.ones((T, 1))

    c = ((1 / (alpha * T)) ** kappa - (1 / (alpha * T)) ** (-kappa)) / (2 * kappa)

    constraints = [
        cp.sum(Z) == 1,
        cp.sum(nu - tau) / (2 * kappa) <= c,
        cp.PowCone3D(nu, ones, Z, 1 / (1 + kappa)),
        cp.PowCone3D(Z, ones, tau, 1 - kappa),
    ]
    risk = Z.T @ (-a)

    objective = cp.Maximize(risk)
    prob = cp.Problem(objective, constraints)

    try:
        if solver in ["CLARABEL", "MOSEK", "SCS"]:
            prob.solve(solver=solver)
        else:
            prob.solve()
    except:
        pass

    if risk.value is None:
        value = None
    else:
        value = risk.value.item()

    if value is None:
        # Primal Formulation
        t = cp.Variable((1, 1))
        z = cp.Variable((1, 1))
        omega = cp.Variable((T, 1))
        psi = cp.Variable((T, 1))
        theta = cp.Variable((T, 1))
        nu = cp.Variable((T, 1))

        ones = np.ones((T, 1))
        constraints = [
            cp.PowCone3D(
                z * (1 + kappa) / (2 * kappa) * ones,
                psi * (1 + kappa) / kappa,
                nu,
                1 / (1 + kappa),
            ),
            cp.PowCone3D(
                omega / (1 - kappa), theta / kappa, -z / (2 * kappa) * ones, (1 - kappa)
            ),
            -a - t + nu + omega <= 0,
            z >= 0,
        ]

        c = ((1 / (alpha * T)) ** kappa - (1 / (alpha * T)) ** (-kappa)) / (2 * kappa)
        risk = t + c * z + cp.sum(psi + theta)

        objective = cp.Minimize(risk)
        prob = cp.Problem(objective, constraints)

        try:
            if solver in ["CLARABEL", "MOSEK", "SCS"]:
                prob.solve(solver=solver)
            else:
                prob.solve(verbose=True)
        except:
            pass

        if risk.value is None:
            value = 0
        else:
            value = risk.value.item()

    return value / 100


def MDD_Abs(X):
    r"""
    Calculate the Maximum Drawdown (MDD) of a returns series
    using uncompounded cumulative returns.

    .. math::
        \text{MDD}(X) = \max_{j \in (0,T)} \left [\max_{t \in (0,j)}
        \left ( \sum_{i=0}^{t}X_{i} \right ) - \sum_{i=0}^{j}X_{i}  \right ]

    Parameters
    ----------
    X : 1d-array
        Returns series, must have Tx1 size.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    Returns
    -------
    value : float
        MDD of an uncompounded cumulative returns.

    """

    a = np.array(X, ndmin=2)
    if a.shape[0] == 1 and a.shape[1] > 1:
        a = a.T
    if a.shape[0] > 1 and a.shape[1] > 1:
        raise ValueError("returns must have Tx1 size")

    prices = np.insert(np.array(a), 0, 1, axis=0)
    NAV = np.cumsum(np.array(prices), axis=0)
    value = 0
    peak = -99999
    for i in NAV:
        if i > peak:
            peak = i
        DD = peak - i
        if DD > value:
            value = DD

    value = np.array(value).item()

    return value


def ADD_Abs(X):
    r"""
    Calculate the Average Drawdown (ADD) of a returns series
    using uncompounded cumulative returns.

    .. math::
        \text{ADD}(X) = \frac{1}{T}\sum_{j=0}^{T}\left [ \max_{t \in (0,j)}
        \left ( \sum_{i=0}^{t}X_{i} \right ) - \sum_{i=0}^{j}X_{i} \right ]

    Parameters
    ----------
    X : 1d-array
        Returns series, must have Tx1 size.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    Returns
    -------
    value : float
        ADD of an uncompounded cumulative returns.

    """

    a = np.array(X, ndmin=2)
    if a.shape[0] == 1 and a.shape[1] > 1:
        a = a.T
    if a.shape[0] > 1 and a.shape[1] > 1:
        raise ValueError("returns must have Tx1 size")

    prices = np.insert(np.array(a), 0, 1, axis=0)
    NAV = np.cumsum(np.array(prices), axis=0)
    value = 0
    peak = -99999
    n = 0
    for i in NAV:
        if i > peak:
            peak = i
        DD = peak - i
        if DD > 0:
            value += DD
        n += 1
    if n == 0:
        value = 0
    else:
        value = value / (n - 1)

    value = np.array(value).item()

    return value


def DaR_Abs(X, alpha=0.05):
    r"""
    Calculate the Drawdown at Risk (DaR) of a returns series
    using uncompounded cumulative returns.

    .. math::
        \text{DaR}_{\alpha}(X) & = \max_{j \in (0,T)} \left \{ \text{DD}(X,j)
        \in \mathbb{R}: F_{\text{DD}} \left ( \text{DD}(X,j) \right )< 1-\alpha
        \right \} \\
        \text{DD}(X,j) & = \max_{t \in (0,j)} \left ( \sum_{i=0}^{t}X_{i}
        \right )- \sum_{i=0}^{j}X_{i}

    Parameters
    ----------
    X : 1d-array
        Returns series, must have Tx1 size..
    alpha : float, optional
        Significance level of DaR. The default is 0.05.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    Returns
    -------
    value : float
        DaR of an uncompounded cumulative returns series.

    """

    a = np.array(X, ndmin=2)
    if a.shape[0] == 1 and a.shape[1] > 1:
        a = a.T
    if a.shape[0] > 1 and a.shape[1] > 1:
        raise ValueError("returns must have Tx1 size")

    prices = np.insert(np.array(a), 0, 1, axis=0)
    NAV = np.cumsum(np.array(prices), axis=0)
    DD = []
    peak = -99999
    for i in NAV:
        if i > peak:
            peak = i
        DD.append(-(peak - i))
    del DD[0]
    sorted_DD = np.sort(np.array(DD), axis=0)
    index = int(np.ceil(alpha * len(sorted_DD)) - 1)
    value = -sorted_DD[index]
    value = np.array(value).item()

    return value


def CDaR_Abs(X, alpha=0.05):
    r"""
    Calculate the Conditional Drawdown at Risk (CDaR) of a returns series
    using uncompounded cumulative returns.

    .. math::
        \text{CDaR}_{\alpha}(X) = \text{DaR}_{\alpha}(X) + \frac{1}{\alpha T}
        \sum_{j=0}^{T} \max \left [ \max_{t \in (0,j)}
        \left ( \sum_{i=0}^{t}X_{i} \right ) - \sum_{i=0}^{j}X_{i}
        - \text{DaR}_{\alpha}(X), 0 \right ]

    Where:

    :math:`\text{DaR}_{\alpha}` is the Drawdown at Risk of an uncompounded
    cumulated return series :math:`X`.

    Parameters
    ----------
    X : 1d-array
        Returns series, must have Tx1 size..
    alpha : float, optional
        Significance level of CDaR. The default is 0.05.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    Returns
    -------
    value : float
        CDaR of an uncompounded cumulative returns series.

    """

    a = np.array(X, ndmin=2)
    if a.shape[0] == 1 and a.shape[1] > 1:
        a = a.T
    if a.shape[0] > 1 and a.shape[1] > 1:
        raise ValueError("returns must have Tx1 size")

    prices = np.insert(np.array(a), 0, 1, axis=0)
    NAV = np.cumsum(np.array(prices), axis=0)
    DD = []
    peak = -99999
    for i in NAV:
        if i > peak:
            peak = i
        DD.append(-(peak - i))
    del DD[0]
    sorted_DD = np.sort(np.array(DD), axis=0)
    index = int(np.ceil(alpha * len(sorted_DD)) - 1)
    sum_var = 0
    for i in range(0, index + 1):
        sum_var = sum_var + sorted_DD[i] - sorted_DD[index]
    value = -sorted_DD[index] - sum_var / (alpha * len(sorted_DD))
    value = np.array(value).item()

    return value


def EDaR_Abs(X, alpha=0.05, solver="CLARABEL"):
    r"""
    Calculate the Entropic Drawdown at Risk (EDaR) of a returns series
    using uncompounded cumulative returns.

    .. math::
        \text{EDaR}_{\alpha}(X) & = \inf_{z>0} \left \{ z
        \ln \left (\frac{M_{\text{DD}(X)}(z^{-1})}{\alpha} \right ) \right \}  \\
        \text{DD}(X,j) & = \max_{t \in (0,j)} \left ( \sum_{i=0}^{t}X_{i}
        \right )- \sum_{i=0}^{j}X_{i} \\

    Parameters
    ----------
    X : 1d-array
        Returns series, must have Tx1 size..
    alpha : float, optional
        Significance level of EDaR. The default is 0.05.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    Returns
    -------
    (value, z) : tuple
        EDaR of an uncompounded cumulative returns series
        and value of z that minimize EDaR.

    """

    a = np.array(X, ndmin=2)
    if a.shape[0] == 1 and a.shape[1] > 1:
        a = a.T
    if a.shape[0] > 1 and a.shape[1] > 1:
        raise ValueError("returns must have Tx1 size")

    prices = np.insert(np.array(a), 0, 1, axis=0)
    NAV = np.cumsum(np.array(prices), axis=0)
    DD = []
    peak = -99999
    for i in NAV:
        if i > peak:
            peak = i
        DD.append(-(peak - i))
    del DD[0]

    (value, t) = EVaR_Hist(np.array(DD), alpha=alpha, solver=solver)

    return (value, t)


def RLDaR_Abs(X, alpha=0.05, kappa=0.3, solver="CLARABEL"):
    r"""
    Calculate the Relativistic Drawdown at Risk (RLDaR) of a returns series
    using uncompounded cumulative returns. I recommend only use this function with MOSEK solver.

    .. math::
        \text{RLDaR}^{\kappa}_{\alpha}(X) & = \text{RLVaR}^{\kappa}_{\alpha}(\text{DD}(X)) \\
        \text{DD}(X,j) & = \max_{t \in (0,j)} \left ( \sum_{i=0}^{t}X_{i}
        \right )- \sum_{i=0}^{j}X_{i} \\

    Parameters
    ----------
    X : 1d-array
        Returns series, must have Tx1 size.
    alpha : float, optional
        Significance level of EVaR. The default is 0.05.
    kappa : float, optional
        Deformation parameter of RLDaR, must be between 0 and 1. The default is 0.3.
    solver: str, optional
        Solver available for CVXPY that supports power cone programming. Used
        to calculate RLVaR, RVRG and RLDaR. The default value is 'CLARABEL'.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    Returns
    -------
    value : tuple
        RLDaR of an uncompounded cumulative returns series.

    """

    a = np.array(X, ndmin=2)
    if a.shape[0] == 1 and a.shape[1] > 1:
        a = a.T
    if a.shape[0] > 1 and a.shape[1] > 1:
        raise ValueError("returns must have Tx1 size")

    prices = np.insert(np.array(a), 0, 1, axis=0)
    NAV = np.cumsum(np.array(prices), axis=0)
    DD = []
    peak = -99999
    for i in NAV:
        if i > peak:
            peak = i
        DD.append(-(peak - i))
    del DD[0]

    value = RLVaR_Hist(np.array(DD), alpha=alpha, kappa=kappa, solver=solver)

    return value


def UCI_Abs(X):
    r"""
    Calculate the Ulcer Index (UCI) of a returns series
    using uncompounded cumulative returns.

    .. math::
        \text{UCI}(X) =\sqrt{\frac{1}{T}\sum_{j=0}^{T} \left [ \max_{t \in
        (0,j)} \left ( \sum_{i=0}^{t}X_{i} \right ) - \sum_{i=0}^{j}X_{i}
        \right ] ^2}

    Parameters
    ----------
    X : 1d-array
        Returns series, must have Tx1 size.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    Returns
    -------
    value : float
        Ulcer Index of an uncompounded cumulative returns.

    """

    a = np.array(X, ndmin=2)
    if a.shape[0] == 1 and a.shape[1] > 1:
        a = a.T
    if a.shape[0] > 1 and a.shape[1] > 1:
        raise ValueError("returns must have Tx1 size")

    prices = np.insert(np.array(a), 0, 1, axis=0)
    NAV = np.cumsum(np.array(prices), axis=0)
    value = 0
    peak = -99999
    n = 0
    for i in NAV:
        if i > peak:
            peak = i
        DD = peak - i
        if DD > 0:
            value += DD**2
        n += 1
    if n == 0:
        value = 0
    else:
        value = np.sqrt(value / (n - 1))

    value = np.array(value).item()

    return value


def MDD_Rel(X):
    r"""
    Calculate the Maximum Drawdown (MDD) of a returns series
    using cumpounded cumulative returns.

    .. math::
        \text{MDD}(X) = \max_{j \in (0,T)}\left[\max_{t \in (0,j)}
        \left ( \prod_{i=0}^{t}(1+X_{i}) \right ) - \prod_{i=0}^{j}(1+X_{i})
        \right]

    Parameters
    ----------
    X : 1d-array
        Returns series, must have Tx1 size.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    Returns
    -------
    value : float
        MDD of a cumpounded cumulative returns.

    """

    a = np.array(X, ndmin=2)
    if a.shape[0] == 1 and a.shape[1] > 1:
        a = a.T
    if a.shape[0] > 1 and a.shape[1] > 1:
        raise ValueError("returns must have Tx1 size")

    prices = 1 + np.insert(np.array(a), 0, 0, axis=0)
    NAV = np.cumprod(prices, axis=0)
    value = 0
    peak = -99999
    for i in NAV:
        if i > peak:
            peak = i
        DD = (peak - i) / peak
        if DD > value:
            value = DD

    value = np.array(value).item()

    return value


def ADD_Rel(X):
    r"""
    Calculate the Average Drawdown (ADD) of a returns series
    using cumpounded cumulative returns.

    .. math::
        \text{ADD}(X) = \frac{1}{T}\sum_{j=0}^{T} \left [ \max_{t \in (0,j)}
        \left ( \prod_{i=0}^{t}(1+X_{i}) \right )- \prod_{i=0}^{j}(1+X_{i})
        \right ]

    Parameters
    ----------
    X : 1d-array
        Returns series, must have Tx1 size.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    Returns
    -------
    value : float
        ADD of a cumpounded cumulative returns.

    """

    a = np.array(X, ndmin=2)
    if a.shape[0] == 1 and a.shape[1] > 1:
        a = a.T
    if a.shape[0] > 1 and a.shape[1] > 1:
        raise ValueError("returns must have Tx1 size")

    prices = 1 + np.insert(np.array(a), 0, 0, axis=0)
    NAV = np.cumprod(prices, axis=0)
    value = 0
    peak = -99999
    n = 0
    for i in NAV:
        if i > peak:
            peak = i
        DD = (peak - i) / peak
        if DD > 0:
            value += DD
        n += 1
    if n == 0:
        value = 0
    else:
        value = value / (n - 1)

    value = np.array(value).item()

    return value


def DaR_Rel(X, alpha=0.05):
    r"""
    Calculate the Drawdown at Risk (DaR) of a returns series
    using cumpounded cumulative returns.

    .. math::
        \text{DaR}_{\alpha}(X) & = \max_{j \in (0,T)} \left \{ \text{DD}(X,j)
        \in \mathbb{R}: F_{\text{DD}} \left ( \text{DD}(X,j) \right )< 1 - \alpha
        \right \} \\
        \text{DD}(X,j) & = \max_{t \in (0,j)} \left ( \prod_{i=0}^{t}(1+X_{i})
        \right )- \prod_{i=0}^{j}(1+X_{i})

    Parameters
    ----------
    X : 1d-array
        Returns series, must have Tx1 size..
    alpha : float, optional
        Significance level of DaR. The default is 0.05.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    Returns
    -------
    value : float
        DaR of a cumpounded cumulative returns series.

    """

    a = np.array(X, ndmin=2)
    if a.shape[0] == 1 and a.shape[1] > 1:
        a = a.T
    if a.shape[0] > 1 and a.shape[1] > 1:
        raise ValueError("X must have Tx1 size")

    prices = 1 + np.insert(np.array(a), 0, 0, axis=0)
    NAV = np.cumprod(prices, axis=0)
    DD = []
    peak = -99999
    for i in NAV:
        if i > peak:
            peak = i
        DD.append(-(peak - i) / peak)
    del DD[0]
    sorted_DD = np.sort(np.array(DD), axis=0)
    index = int(np.ceil(alpha * len(sorted_DD)) - 1)
    value = -sorted_DD[index]
    value = np.array(value).item()

    return value


def CDaR_Rel(X, alpha=0.05):
    r"""
    Calculate the Conditional Drawdown at Risk (CDaR) of a returns series
    using cumpounded cumulative returns.

    .. math::
        \text{CDaR}_{\alpha}(X) = \text{DaR}_{\alpha}(X) + \frac{1}{\alpha T}
        \sum_{i=0}^{T} \max \left [ \max_{t \in (0,T)}
        \left ( \prod_{i=0}^{t}(1+X_{i}) \right )- \prod_{i=0}^{j}(1+X_{i})
        - \text{DaR}_{\alpha}(X), 0 \right ]

    Where:

    :math:`\text{DaR}_{\alpha}` is the Drawdown at Risk of a cumpound
    cumulated return series :math:`X`.

    Parameters
    ----------
    X : 1d-array
        Returns series, must have Tx1 size..
    alpha : float, optional
        Significance level of CDaR. The default is 0.05.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    Returns
    -------
    value : float
        CDaR of a cumpounded cumulative returns series.

    """

    a = np.array(X, ndmin=2)
    if a.shape[0] == 1 and a.shape[1] > 1:
        a = a.T
    if a.shape[0] > 1 and a.shape[1] > 1:
        raise ValueError("X must have Tx1 size")

    prices = 1 + np.insert(np.array(a), 0, 0, axis=0)
    NAV = np.cumprod(prices, axis=0)
    DD = []
    peak = -99999
    for i in NAV:
        if i > peak:
            peak = i
        DD.append(-(peak - i) / peak)
    del DD[0]
    sorted_DD = np.sort(np.array(DD), axis=0)
    index = int(np.ceil(alpha * len(sorted_DD)) - 1)
    sum_var = 0
    for i in range(0, index + 1):
        sum_var = sum_var + sorted_DD[i] - sorted_DD[index]
    value = -sorted_DD[index] - sum_var / (alpha * len(sorted_DD))
    value = np.array(value).item()

    return value


def EDaR_Rel(X, alpha=0.05, solver="CLARABEL"):
    r"""
    Calculate the Entropic Drawdown at Risk (EDaR) of a returns series
    using cumpounded cumulative returns.

    .. math::
        \text{EDaR}_{\alpha}(X) & = \inf_{z>0} \left \{ z
        \ln \left (\frac{M_{\text{DD}(X)}(z^{-1})}{\alpha} \right ) \right \}  \\
        \text{DD}(X,j) & = \max_{t \in (0,j)} \left ( \prod_{i=0}^{t}(1+X_{i})
        \right )- \prod_{i=0}^{j}(1+X_{i})

    Parameters
    ----------
    X : 1d-array
        Returns series, must have Tx1 size..
    alpha : float, optional
        Significance level of EDaR. The default is 0.05.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    Returns
    -------
    (value, z) : tuple
        EDaR of a cumpounded cumulative returns series
        and value of z that minimize EDaR.

    """

    a = np.array(X, ndmin=2)
    if a.shape[0] == 1 and a.shape[1] > 1:
        a = a.T
    if a.shape[0] > 1 and a.shape[1] > 1:
        raise ValueError("X must have Tx1 size")

    prices = 1 + np.insert(np.array(a), 0, 0, axis=0)
    NAV = np.cumprod(prices, axis=0)
    DD = []
    peak = -99999
    for i in NAV:
        if i > peak:
            peak = i
        DD.append(-(peak - i) / peak)
    del DD[0]

    (value, t) = EVaR_Hist(np.array(DD), alpha=alpha, solver=solver)

    return (value, t)


def RLDaR_Rel(X, alpha=0.05, kappa=0.3, solver="CLARABEL"):
    r"""
    Calculate the Relativistic Drawdown at Risk (RLDaR) of a returns series
    using compounded cumulative returns. I recommend only use this function with MOSEK solver.

    .. math::
        \text{RLDaR}^{\kappa}_{\alpha}(X) & = \text{RLVaR}^{\kappa}_{\alpha}(\text{DD}(X)) \\
        \text{DD}(X,j) & = \max_{t \in (0,j)} \left ( \prod_{i=0}^{t}(1+X_{i})
        \right )- \prod_{i=0}^{j}(1+X_{i}) \\

    Parameters
    ----------
    X : 1d-array
        Returns series, must have Tx1 size.
    alpha : float, optional
        Significance level of RLDaR. The default is 0.05.
    kappa : float, optional
        Deformation parameter of RLDaR, must be between 0 and 1. The default is 0.3.
    solver: str, optional
        Solver available for CVXPY that supports power cone programming. Used
        to calculate RLVaR, RVRG and RLDaR. The default value is 'CLARABEL'.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    Returns
    -------
    value : tuple
        RLDaR of a compounded cumulative returns series.

    """

    a = np.array(X, ndmin=2)
    if a.shape[0] == 1 and a.shape[1] > 1:
        a = a.T
    if a.shape[0] > 1 and a.shape[1] > 1:
        raise ValueError("X must have Tx1 size")

    prices = 1 + np.insert(np.array(a), 0, 0, axis=0)
    NAV = np.cumprod(prices, axis=0)
    DD = []
    peak = -99999
    for i in NAV:
        if i > peak:
            peak = i
        DD.append(-(peak - i) / peak)
    del DD[0]

    value = RLVaR_Hist(np.array(DD), alpha=alpha, kappa=kappa, solver=solver)

    return value


def UCI_Rel(X):
    r"""
    Calculate the Ulcer Index (UCI) of a returns series
    using cumpounded cumulative returns.

    .. math::
        \text{UCI}(X) =\sqrt{\frac{1}{T}\sum_{j=0}^{T} \left [ \max_{t \in
        (0,j)} \left ( \prod_{i=0}^{t}(1+X_{i}) \right )- \prod_{i=0}^{j}
        (1+X_{i}) \right ] ^2}

    Parameters
    ----------
    X : 1d-array
        Returns series, must have Tx1 size.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    Returns
    -------
    value : float
        Ulcer Index of a cumpounded cumulative returns.

    """

    a = np.array(X, ndmin=2)
    if a.shape[0] == 1 and a.shape[1] > 1:
        a = a.T
    if a.shape[0] > 1 and a.shape[1] > 1:
        raise ValueError("returns must have Tx1 size")

    prices = 1 + np.insert(np.array(a), 0, 0, axis=0)
    NAV = np.cumprod(prices, axis=0)
    value = 0
    peak = -99999
    n = 0
    for i in NAV:
        if i > peak:
            peak = i
        DD = (peak - i) / peak
        if DD > 0:
            value += DD**2
        n += 1
    if n == 0:
        value = 0
    else:
        value = np.sqrt(value / (n - 1))

    value = np.array(value).item()

    return value


def GMD(X):
    r"""
    Calculate the Gini Mean Difference (GMD) of a returns series.

    Parameters
    ----------
    X : 1d-array
        Returns series, must have Tx1 size.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    Returns
    -------
    value : float
        Gini Mean Difference of a returns series.

    """

    a = np.array(X, ndmin=2)
    if a.shape[0] == 1 and a.shape[1] > 1:
        a = a.T
    if a.shape[0] > 1 and a.shape[1] > 1:
        raise ValueError("returns must have Tx1 size")

    T = a.shape[0]
    w_ = owa.owa_gmd(T)
    value = (w_.T @ np.sort(a, axis=0)).item()

    return value


def TG(X, alpha=0.05, a_sim=100):
    r"""
    Calculate the Tail Gini of a returns series.

    Parameters
    ----------
    X : 1d-array
        Returns series, must have Tx1 size.
    alpha : float, optional
        Significance level of Tail Gini. The default is 0.05.
    a_sim : float, optional
        Number of CVaRs used to approximate Tail Gini. The default is 100.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    Returns
    -------
    value : float
        Ulcer Index of a cumpounded cumulative returns.

    """

    a = np.array(X, ndmin=2)
    if a.shape[0] == 1 and a.shape[1] > 1:
        a = a.T
    if a.shape[0] > 1 and a.shape[1] > 1:
        raise ValueError("returns must have Tx1 size")

    T = a.shape[0]
    w_ = owa.owa_tg(T, alpha, a_sim)
    value = (w_.T @ np.sort(a, axis=0)).item()

    return value


def RG(X):
    r"""
    Calculate the range of a returns series.

    Parameters
    ----------
    X : 1d-array
        Returns series, must have Tx1 size.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    Returns
    -------
    value : float
        Ulcer Index of a cumpounded cumulative returns.

    """

    a = np.array(X, ndmin=2)
    if a.shape[0] == 1 and a.shape[1] > 1:
        a = a.T
    if a.shape[0] > 1 and a.shape[1] > 1:
        raise ValueError("returns must have Tx1 size")

    T = a.shape[0]
    w_ = owa.owa_rg(T)
    value = (w_.T @ np.sort(a, axis=0)).item()

    return value


def VRG(X, alpha=0.05, beta=None):
    r"""
    Calculate the CVaR range of a returns series.

    Parameters
    ----------
    X : 1d-array
        Returns series, must have Tx1 size.
    alpha : float, optional
        Significance level of VaR of losses. The default is 0.05.
    beta : float, optional
        Significance level of VaR of gains. If None it duplicates alpha value.
        The default is None.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    Returns
    -------
    value : float
        Ulcer Index of a cumpounded cumulative returns.

    """

    a = np.array(X, ndmin=2)
    if a.shape[0] == 1 and a.shape[1] > 1:
        a = a.T
    if a.shape[0] > 1 and a.shape[1] > 1:
        raise ValueError("returns must have Tx1 size")

    if beta is None:
        beta = alpha

    value_L = VaR_Hist(a, alpha=alpha)
    value_G = VaR_Hist(-a, alpha=beta)

    value = value_L + value_G

    return value


def CVRG(X, alpha=0.05, beta=None):
    r"""
    Calculate the CVaR range of a returns series.

    Parameters
    ----------
    X : 1d-array
        Returns series, must have Tx1 size.
    alpha : float, optional
        Significance level of CVaR of losses. The default is 0.05.
    beta : float, optional
        Significance level of CVaR of gains. If None it duplicates alpha value.
        The default is None.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    Returns
    -------
    value : float
        Ulcer Index of a cumpounded cumulative returns.

    """

    a = np.array(X, ndmin=2)
    if a.shape[0] == 1 and a.shape[1] > 1:
        a = a.T
    if a.shape[0] > 1 and a.shape[1] > 1:
        raise ValueError("returns must have Tx1 size")

    T = a.shape[0]
    w_ = owa.owa_cvrg(T, alpha=alpha, beta=beta)
    value = (w_.T @ np.sort(a, axis=0)).item()

    return value


def TGRG(X, alpha=0.05, a_sim=100, beta=None, b_sim=None):
    r"""
    Calculate the Tail Gini range of a returns series.

    Parameters
    ----------
    X : 1d-array
        Returns series, must have Tx1 size.
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

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    Returns
    -------
    value : float
        Ulcer Index of a cumpounded cumulative returns.

    """

    a = np.array(X, ndmin=2)
    if a.shape[0] == 1 and a.shape[1] > 1:
        a = a.T
    if a.shape[0] > 1 and a.shape[1] > 1:
        raise ValueError("returns must have Tx1 size")

    T = a.shape[0]
    w_ = owa.owa_tgrg(T, alpha=alpha, a_sim=a_sim, beta=beta, b_sim=b_sim)
    value = (w_.T @ np.sort(a, axis=0)).item()

    return value


def EVRG(X, alpha=0.05, beta=None, solver="CLARABEL"):
    r"""
    Calculate the CVaR range of a returns series.

    Parameters
    ----------
    X : 1d-array
        Returns series, must have Tx1 size.
    alpha : float, optional
        Significance level of EVaR of losses. The default is 0.05.
    beta : float, optional
        Significance level of EVaR of gains. If None it duplicates alpha value.
        The default is None.
    solver: str, optional
        Solver available for CVXPY that supports exponential cone programming.
        Used to calculate EVaR, EVRG and EDaR. The default value is 'CLARABEL'.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    Returns
    -------
    value : float
        Ulcer Index of a cumpounded cumulative returns.

    """

    a = np.array(X, ndmin=2)
    if a.shape[0] == 1 and a.shape[1] > 1:
        a = a.T
    if a.shape[0] > 1 and a.shape[1] > 1:
        raise ValueError("returns must have Tx1 size")

    if beta is None:
        beta = alpha

    value_L = EVaR_Hist(a, alpha=alpha, solver=solver)[0]
    value_G = EVaR_Hist(-a, alpha=beta, solver=solver)[0]

    value = value_L + value_G

    return value


def RVRG(X, alpha=0.05, beta=None, kappa=0.3, kappa_g=None, solver="CLARABEL"):
    r"""
    Calculate the CVaR range of a returns series.

    Parameters
    ----------
    X : 1d-array
        Returns series, must have Tx1 size.
    alpha : float, optional
        Significance level of RLVaR of losses. The default is 0.05.
    beta : float, optional
        Significance level of RLVaR of gains. If None it duplicates alpha value.
        The default is None.
    kappa : float, optional
        Deformation parameter of RLVaR for losses, must be between 0 and 1.
        The default is 0.3.
    kappa_g : float, optional
        Deformation parameter of RLVaR for gains, must be between 0 and 1.
        The default is None.
    solver: str, optional
        Solver available for CVXPY that supports power cone programming.
        Used to calculate EVaR, EVRG and EDaR. The default value is 'CLARABEL'.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    Returns
    -------
    value : float
        Ulcer Index of a cumpounded cumulative returns.

    """

    a = np.array(X, ndmin=2)
    if a.shape[0] == 1 and a.shape[1] > 1:
        a = a.T
    if a.shape[0] > 1 and a.shape[1] > 1:
        raise ValueError("returns must have Tx1 size")

    if beta is None:
        beta = alpha
    if kappa_g is None:
        kappa_g = kappa

    value_L = RLVaR_Hist(a, alpha=alpha, kappa=kappa, solver=solver)
    value_G = RLVaR_Hist(-a, alpha=beta, kappa=kappa_g, solver=solver)

    value = value_L + value_G

    return value


def L_Moment(X, k=2):
    r"""
    Calculate the kth l-moment of a returns series.

    .. math:
        \lambda_k = {\tbinom{T}{k}}^{-1} \mathop{\sum \sum \ldots \sum}_{1
        \leq i_{1} < i_{2} \cdots < i_{k} \leq n} \frac{1}{k}
        \sum^{k-1}_{j=0} (-1)^{j} \binom{k-1}{j} y_{[i_{k-j}]} \\

    Where $y_{[i]}$ is the ith-ordered statistic.

    Parameters
    ----------
    X : 1d-array
        Returns series, must have Tx1 size.
    k : int
        Order of the l-moment. Must be an integer higher or equal than 1.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    Returns
    -------
    value : float
        Kth l-moment of a returns series.

    """

    a = np.array(X, ndmin=2)
    if a.shape[0] == 1 and a.shape[1] > 1:
        a = a.T
    if a.shape[0] > 1 and a.shape[1] > 1:
        raise ValueError("returns must have Tx1 size")

    T = a.shape[0]
    w_ = owa.owa_l_moment(T, k=k)
    value = (w_.T @ np.sort(a, axis=0)).item()

    return value


def L_Moment_CRM(X, k=4, method="MSD", g=0.5, max_phi=0.5, solver="CLARABEL"):
    r"""
    Calculate a custom convex risk measure that is a weighted average of
    first k-th l-moments.

    Parameters
    ----------
    X : 1d-array
        Returns series, must have Tx1 size.
    k : int
        Order of the l-moment. Must be an integer higher or equal than 2.
    method : str, optional
        Method to calculate the weights used to combine the l-moments with
        order higher than 2. The default value is 'MSD'. Possible values are:

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
        The default value is None.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    Returns
    -------
    value : float
        Custom convex risk measure that is a weighted average of first k-th l-moments of a returns series.

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

    a = np.array(X, ndmin=2)
    if a.shape[0] == 1 and a.shape[1] > 1:
        a = a.T
    if a.shape[0] > 1 and a.shape[1] > 1:
        raise ValueError("returns must have Tx1 size")

    T = a.shape[0]
    w_ = owa.owa_l_moment_crm(
        T, k=k, method=method, g=g, max_phi=max_phi, solver=solver
    )
    value = (w_.T @ np.sort(a, axis=0)).item()

    return value


def NEA(w):
    r"""
    Calculate the number of effective assets (NEA) that is the inverse of the
    Herfindahl Hirschman index (HHI).

    Parameters
    ----------
    w : DataFrame or Series of shape (n_assets, 1)
        Portfolio weights, where n_assets is the number of assets.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    Returns
    -------
    value : float
        The NEA of the portfolio.
    """

    a = np.array(w, ndmin=2)
    if a.shape[0] == 1 and a.shape[1] > 1:
        a = a.T
    if a.shape[0] > 1 and a.shape[1] > 1:
        raise ValueError("w must have n_assets x 1 size")

    value = 1 / np.sum(a**2)

    return value


###############################################################################
# Risk Adjusted Return Ratios
###############################################################################


def Sharpe_Risk(
    returns,
    w=None,
    cov=None,
    rm="MV",
    rf=0,
    alpha=0.05,
    a_sim=100,
    beta=None,
    b_sim=None,
    kappa=0.3,
    kappa_g=None,
    solver="CLARABEL",
):
    r"""
    Calculate the risk measure available on the Sharpe function.

    Parameters
    ----------
    w : DataFrame or 1d-array of shape (n_assets, 1)
        Weights matrix, where n_assets is the number of assets.
    cov : DataFrame of shape (n_assets, n_assets)
        Covariance matrix, where n_assets is the number of assets.
    returns : DataFrame or nd-array of shape (n_samples, n_features)
        Features matrix, where n_samples is the number of samples and
        n_features is the number of features.
    rm : str, optional
        Risk measure used in the denominator of the ratio. The default is
        'MV'. Possible values are:

        - 'MV': Standard Deviation.
        - 'KT': Square Root Kurtosis.
        - 'MAD': Mean Absolute Deviation.
        - 'GMD': Gini Mean Difference.
        - 'MSV': Semi Standard Deviation.
        - 'SKT': Square Root Semi Kurtosis.
        - 'FLPM': First Lower Partial Moment (Omega Ratio).
        - 'SLPM': Second Lower Partial Moment (Sortino Ratio).
        - 'VaR': Value at Risk.
        - 'CVaR': Conditional Value at Risk.
        - 'TG': Tail Gini.
        - 'EVaR': Entropic Value at Risk.
        - 'RLVaR': Relativistic Value at Risk. I recommend only use this function with MOSEK solver.
        - 'WR': Worst Realization (Minimax).
        - 'RG': Range of returns.
        - 'VRG' VaR range of returns.
        - 'CVRG': CVaR range of returns.
        - 'TGRG': Tail Gini range of returns.
        - 'EVRG': EVaR range of returns.
        - 'RVRG': RLVaR range of returns. I recommend only use this function with MOSEK solver.
        - 'MDD': Maximum Drawdown of uncompounded cumulative returns (Calmar Ratio).
        - 'ADD': Average Drawdown of uncompounded cumulative returns.
        - 'DaR': Drawdown at Risk of uncompounded cumulative returns.
        - 'CDaR': Conditional Drawdown at Risk of uncompounded cumulative returns.
        - 'EDaR': Entropic Drawdown at Risk of uncompounded cumulative returns.
        - 'RLDaR': Relativistic Drawdown at Risk of uncompounded cumulative returns. I recommend only use this risk measure with MOSEK solver.
        - 'UCI': Ulcer Index of uncompounded cumulative returns.
        - 'MDD_Rel': Maximum Drawdown of compounded cumulative returns (Calmar Ratio).
        - 'ADD_Rel': Average Drawdown of compounded cumulative returns.
        - 'DaR_Rel': Drawdown at Risk of compounded cumulative returns.
        - 'CDaR_Rel': Conditional Drawdown at Risk of compounded cumulative returns.
        - 'EDaR_Rel': Entropic Drawdown at Risk of compounded cumulative returns.
        - 'RLDaR_Rel': Relativistic Drawdown at Risk of compounded cumulative returns. I recommend only use this risk measure with MOSEK solver.
        - 'UCI_Rel': Ulcer Index of compounded cumulative returns.

    rf : float, optional
        Risk free rate. The default is 0.
    alpha : float, optional
        Significance level of VaR, CVaR, EVaR, RLVaR, DaR, CDaR, EDaR, RLDaR
        and Tail Gini of losses. The default is 0.05.
    a_sim : float, optional
        Number of CVaRs used to approximate Tail Gini of losses. The default is 100.
    beta : float, optional
        Significance level of CVaR and Tail Gini of gains. If None it
        duplicates alpha value. The default is None.
    b_sim : float, optional
        Number of CVaRs used to approximate Tail Gini of gains. If None it
        duplicates a_sim value. The default is None.
    kappa : float, optional
        Deformation parameter of RLVaR and RLDaR for losses, must be between 0 and 1.
        The default is 0.3.
    kappa_g : float, optional
        Deformation parameter of RLVaR and RLDaR for gains, must be between 0 and 1.
        The default is None.
    solver: str, optional
        Solver available for CVXPY that supports exponential and power cone
        programming. Used to calculate RLVaR and RLDaR. The default value is
        'CLARABEL'.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    Returns
    -------
    value : float
        Risk measure of the portfolio.

    """

    if isinstance(returns, pd.Series):
        returns_ = returns.to_frame()
    elif isinstance(returns, pd.DataFrame):
        returns_ = returns.to_numpy()
    else:
        returns_ = np.array(returns, ndmin=2)

    if returns_.shape[1] == 1:
        w_ = np.array([[1]])
    else:
        if w is None:
            raise ValueError("weights must have n_assets x 1 size")
        else:
            w_ = np.array(w, ndmin=2)

    if w_.shape[0] == 1 and w_.shape[1] > 1:
        w_ = w_.T
    if w_.shape[0] > 1 and w_.shape[1] > 1:
        raise ValueError("weights must have n_assets x 1 size")

    if cov is None:
        cov_ = np.array(np.cov(returns_, rowvar=False), ndmin=2)
    else:
        cov_ = np.array(cov, ndmin=2)

    a = returns_ @ w_
    if rm == "MV":
        risk = w_.T @ cov_ @ w_
        risk = np.sqrt(risk.item())
    elif rm == "MAD":
        risk = MAD(a)
    elif rm == "GMD":
        risk = GMD(a)
    elif rm == "MSV":
        risk = SemiDeviation(a)
    elif rm == "FLPM":
        risk = LPM(a, MAR=rf, p=1)
    elif rm == "SLPM":
        risk = LPM(a, MAR=rf, p=2)
    elif rm == "VaR":
        risk = VaR_Hist(a, alpha=alpha)
    elif rm == "CVaR":
        risk = CVaR_Hist(a, alpha=alpha)
    elif rm == "TG":
        risk = TG(a, alpha=alpha, a_sim=a_sim)
    elif rm == "EVaR":
        risk = EVaR_Hist(a, alpha=alpha, solver=solver)[0]
    elif rm == "RLVaR":
        risk = RLVaR_Hist(a, alpha=alpha, kappa=kappa, solver=solver)
    elif rm == "WR":
        risk = WR(a)
    elif rm == "RG":
        risk = RG(a)
    elif rm == "VRG":
        risk = VRG(a, alpha=alpha, beta=beta)
    elif rm == "CVRG":
        risk = CVRG(a, alpha=alpha, beta=beta)
    elif rm == "TGRG":
        risk = TGRG(a, alpha=alpha, a_sim=a_sim, beta=beta, b_sim=b_sim)
    elif rm == "EVRG":
        risk = EVRG(a, alpha=alpha, beta=beta, solver=solver)
    elif rm == "RVRG":
        risk = RVRG(
            a, alpha=alpha, beta=beta, kappa=kappa, kappa_g=kappa_g, solver=solver
        )
    elif rm == "MDD":
        risk = MDD_Abs(a)
    elif rm == "ADD":
        risk = ADD_Abs(a)
    elif rm == "DaR":
        risk = DaR_Abs(a, alpha=alpha)
    elif rm == "CDaR":
        risk = CDaR_Abs(a, alpha=alpha)
    elif rm == "EDaR":
        risk = EDaR_Abs(a, alpha=alpha)[0]
    elif rm == "RLDaR":
        risk = RLDaR_Abs(a, alpha=alpha, kappa=kappa, solver=solver)
    elif rm == "UCI":
        risk = UCI_Abs(a)
    elif rm == "MDD_Rel":
        risk = MDD_Rel(a)
    elif rm == "ADD_Rel":
        risk = ADD_Rel(a)
    elif rm == "DaR_Rel":
        risk = DaR_Rel(a, alpha=alpha)
    elif rm == "CDaR_Rel":
        risk = CDaR_Rel(a, alpha=alpha)
    elif rm == "EDaR_Rel":
        risk = EDaR_Rel(a, alpha=alpha)[0]
    elif rm == "RLDaR_Rel":
        risk = RLDaR_Rel(a, alpha=alpha, kappa=kappa, solver=solver)
    elif rm == "UCI_Rel":
        risk = UCI_Rel(a)
    elif rm == "KT":
        risk = Kurtosis(a)
    elif rm == "SKT":
        risk = SemiKurtosis(a)

    value = risk

    return value


def Sharpe(
    returns,
    w=None,
    mu=None,
    cov=None,
    rm="MV",
    rf=0,
    alpha=0.05,
    a_sim=100,
    beta=None,
    b_sim=None,
    kappa=0.3,
    kappa_g=None,
    solver="CLARABEL",
):
    r"""
    Calculate the Risk Adjusted Return Ratio from a portfolio returns series.

    .. math::
        \text{Sharpe}(X) =  \frac{\mathbb{E}(X) -
        r_{f}}{\phi(X)}

    Where:

    :math:`X` is the vector of portfolio returns.

    :math:`r_{f}` is the risk free rate, when the risk measure is

    :math:`\text{LPM}` uses instead of :math:`r_{f}` the :math:`\text{MAR}`.

    :math:`\phi(X)` is a convex risk measure. The risk measures availabe are:

    Parameters
    ----------

    returns : DataFrame or nd-array of shape (n_samples, n_features)
        Features matrix, where n_samples is the number of samples and
        n_features is the number of features.
    w : DataFrame or 1d-array of shape (n_assets, 1)
        Weights matrix, where n_assets is the number of assets.
    mu : DataFrame or nd-array of shape (1, n_assets)
        Vector of expected returns, where n_assets is the number of assets.
    cov : DataFrame of shape (n_assets, n_assets)
        Covariance matrix, where n_assets is the number of assets.
    rm : str, optional
        Risk measure used in the denominator of the ratio. The default is
        'MV'. Possible values are:

        - 'MV': Standard Deviation.
        - 'KT': Square Root Kurtosis.
        - 'MAD': Mean Absolute Deviation.
        - 'GMD': Gini Mean Difference.
        - 'MSV': Semi Standard Deviation.
        - 'SKT': Square Root Semi Kurtosis.
        - 'FLPM': First Lower Partial Moment (Omega Ratio).
        - 'SLPM': Second Lower Partial Moment (Sortino Ratio).
        - 'VaR': Value at Risk.
        - 'CVaR': Conditional Value at Risk.
        - 'TG': Tail Gini.
        - 'EVaR': Entropic Value at Risk.
        - 'RLVaR': Relativistic Value at Risk. I recommend only use this function with MOSEK solver.
        - 'WR': Worst Realization (Minimax).
        - 'RG': Range of returns.
        - 'VRG' VaR range of returns.
        - 'CVRG': CVaR range of returns.
        - 'TGRG': Tail Gini range of returns.
        - 'EVRG': EVaR range of returns.
        - 'RVRG': RLVaR range of returns. I recommend only use this function with MOSEK solver.
        - 'MDD': Maximum Drawdown of uncompounded cumulative returns (Calmar Ratio).
        - 'ADD': Average Drawdown of uncompounded cumulative returns.
        - 'DaR': Drawdown at Risk of uncompounded cumulative returns.
        - 'CDaR': Conditional Drawdown at Risk of uncompounded cumulative returns.
        - 'EDaR': Entropic Drawdown at Risk of uncompounded cumulative returns.
        - 'RLDaR': Relativistic Drawdown at Risk of uncompounded cumulative returns. I recommend only use this function with MOSEK solver.
        - 'UCI': Ulcer Index of uncompounded cumulative returns.
        - 'MDD_Rel': Maximum Drawdown of compounded cumulative returns (Calmar Ratio).
        - 'ADD_Rel': Average Drawdown of compounded cumulative returns.
        - 'DaR_Rel': Drawdown at Risk of compounded cumulative returns.
        - 'CDaR_Rel': Conditional Drawdown at Risk of compounded cumulative returns.
        - 'EDaR_Rel': Entropic Drawdown at Risk of compounded cumulative returns.
        - 'RLDaR_Rel': Relativistic Drawdown at Risk of compounded cumulative returns. I recommend only use this function with MOSEK solver.
        - 'UCI_Rel': Ulcer Index of compounded cumulative returns.

    rf : float, optional
        Risk free rate. The default is 0.
    alpha : float, optional
        Significance level of VaR, CVaR, EVaR, RLVaR, DaR, CDaR, EDaR, RLDaR and Tail Gini of losses.
        The default is 0.05.
    a_sim : float, optional
        Number of CVaRs used to approximate Tail Gini of losses. The default is 100.
    beta : float, optional
        Significance level of CVaR and Tail Gini of gains. If None it duplicates alpha value.
        The default is None.
    b_sim : float, optional
        Number of CVaRs used to approximate Tail Gini of gains. If None it duplicates a_sim value.
        The default is None.
    kappa : float, optional
        Deformation parameter of RLVaR and RLDaR for losses, must be between 0 and 1.
        The default is 0.3.
    kappa_g : float, optional
        Deformation parameter of RLVaR and RLDaR for gains, must be between 0 and 1.
        The default is None.
    solver: str, optional
        Solver available for CVXPY that supports power cone programming. Used to calculate RLVaR and RLDaR.
        The default value is None.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    Returns
    -------
    value : float
        Risk adjusted return ratio of :math:`X`.

    """

    if isinstance(returns, pd.Series):
        returns_ = returns.to_frame()
    elif isinstance(returns, pd.DataFrame):
        returns_ = returns.to_numpy()
    else:
        returns_ = np.array(returns, ndmin=2)

    if returns_.shape[1] == 1:
        w_ = np.array([[1]])
    else:
        if w is None:
            raise ValueError("weights must have n_assets x 1 size")
        else:
            w_ = np.array(w, ndmin=2)

    if w_.shape[0] == 1 and w_.shape[1] > 1:
        w_ = w_.T
    if w_.shape[0] > 1 and w_.shape[1] > 1:
        raise ValueError("weights must have n_assets x 1 size")

    if cov is None:
        cov_ = np.array(np.cov(returns_, rowvar=False), ndmin=2)
    else:
        cov_ = np.array(cov, ndmin=2)

    if mu is None:
        mu_ = np.array(np.mean(returns_, axis=0), ndmin=2)
    else:
        mu_ = np.array(mu, ndmin=2)

    ret = mu_ @ w_
    ret = ret.item()

    risk = Sharpe_Risk(
        returns=returns_,
        w=w_,
        cov=cov_,
        rm=rm,
        rf=rf,
        alpha=alpha,
        a_sim=a_sim,
        beta=beta,
        b_sim=b_sim,
        kappa=kappa,
        kappa_g=kappa_g,
        solver=solver,
    )

    value = (ret - rf) / risk

    return value


###############################################################################
# Risk Contribution Vectors
###############################################################################


def Risk_Contribution(
    w,
    returns,
    cov=None,
    rm="MV",
    rf=0,
    alpha=0.05,
    a_sim=100,
    beta=None,
    b_sim=None,
    kappa=0.3,
    kappa_g=None,
    solver="CLARABEL",
):
    r"""
    Calculate the risk contribution for each asset based on the selected risk measure.

    Parameters
    ----------
    w : DataFrame or Series of shape (n_assets, 1)
        Portfolio weights, where n_assets is the number of assets.
    returns : DataFrame or nd-array of shape (n_samples, n_features)
        Features matrix, where n_samples is the number of samples and
        n_features is the number of features.
    cov : DataFrame of shape (n_assets, n_assets)
        Covariance matrix, where n_assets is the number of assets.
    rm : str, optional
        Risk measure used in the denominator of the ratio. The default is
        'MV'. Possible values are:

        - 'MV': Standard Deviation.
        - 'KT': Square Root Kurtosis.
        - 'MAD': Mean Absolute Deviation.
        - 'GMD': Gini Mean Difference.
        - 'MSV': Semi Standard Deviation.
        - 'SKT': Square Root Semi Kurtosis.
        - 'FLPM': First Lower Partial Moment (Omega Ratio).
        - 'SLPM': Second Lower Partial Moment (Sortino Ratio).
        - 'VaR': Value at Risk.
        - 'CVaR': Conditional Value at Risk.
        - 'TG': Tail Gini.
        - 'EVaR': Entropic Value at Risk.
        - 'RLVaR': Relativistic Value at Risk. I recommend only use this function with MOSEK solver.
        - 'WR': Worst Realization (Minimax).
        - 'RG': Range of returns.
        - 'VRG' VaR range of returns.
        - 'CVRG': CVaR range of returns.
        - 'TGRG': Tail Gini range of returns.
        - 'EVRG': EVaR range of returns.
        - 'RVRG': RLVaR range of returns. I recommend only use this function with MOSEK solver.
        - 'MDD': Maximum Drawdown of uncompounded cumulative returns (Calmar Ratio).
        - 'ADD': Average Drawdown of uncompounded cumulative returns.
        - 'DaR': Drawdown at Risk of uncompounded cumulative returns.
        - 'CDaR': Conditional Drawdown at Risk of uncompounded cumulative returns.
        - 'EDaR': Entropic Drawdown at Risk of uncompounded cumulative returns.
        - 'RLDaR': Relativistic Drawdown at Risk of uncompounded cumulative returns. I recommend only use this function with MOSEK solver.
        - 'UCI': Ulcer Index of uncompounded cumulative returns.
        - 'MDD_Rel': Maximum Drawdown of compounded cumulative returns (Calmar Ratio).
        - 'ADD_Rel': Average Drawdown of compounded cumulative returns.
        - 'CDaR_Rel': Conditional Drawdown at Risk of compounded cumulative returns.
        - 'EDaR_Rel': Entropic Drawdown at Risk of compounded cumulative returns.
        - 'RLDaR_Rel': Relativistic Drawdown at Risk of compounded cumulative returns. I recommend only use this function with MOSEK solver.
        - 'UCI_Rel': Ulcer Index of compounded cumulative returns.

    rf : float, optional
        Risk free rate. The default is 0.
    alpha : float, optional
        Significance level of VaR, CVaR, EVaR, RLVaR, DaR, CDaR, EDaR, RLDaR and Tail Gini of losses.
        The default is 0.05.
    a_sim : float, optional
        Number of CVaRs used to approximate Tail Gini of losses. The default is 100.
    beta : float, optional
        Significance level of CVaR and Tail Gini of gains. If None it duplicates alpha value.
        The default is None.
    b_sim : float, optional
        Number of CVaRs used to approximate Tail Gini of gains. If None it duplicates a_sim value.
        The default is None.
    kappa : float, optional
        Deformation parameter of RLVaR and RLDaR for losses, must be between 0 and 1.
        The default is 0.3.
    kappa_g : float, optional
        Deformation parameter of RLVaR and RLDaR for gains, must be between 0 and 1.
        The default is None.
    solver: str, optional
        Solver available for CVXPY that supports power cone programming. Used to calculate RLVaR and RLDaR.
        The default value is None.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    Returns
    -------
    value : float
        Risk measure of the portfolio.

    """

    w_ = np.array(w, ndmin=2)
    if w_.shape[0] == 1 and w_.shape[1] > 1:
        w_ = w_.T
    if w_.shape[0] > 1 and w_.shape[1] > 1:
        raise ValueError("weights must have n_assets x 1 size")

    if isinstance(returns, pd.Series):
        returns_ = returns.to_frame()
        returns_ = returns.to_numpy()
    elif isinstance(returns, pd.DataFrame):
        returns_ = returns.to_numpy()
    else:
        returns_ = np.array(returns, ndmin=2)

    if cov is None:
        cov_ = np.array(np.cov(returns_, rowvar=False), ndmin=2)
    else:
        cov_ = np.array(cov, ndmin=2)

    RC = []
    if rm in ["EVaR", "EDaR", "RLVaR", "RLDaR", "EVRG", "RVRG"]:
        d_i = 0.0001
    else:
        d_i = 0.0000001

    for i in range(0, w_.shape[0]):
        delta = np.zeros((w_.shape[0], 1))
        delta[i, 0] = d_i
        w_1 = w_ + delta
        w_2 = w_ - delta
        a_1 = returns_ @ w_1
        a_2 = returns_ @ w_2
        if rm == "MV":
            risk_1 = w_1.T @ cov_ @ w_1
            risk_1 = np.sqrt(risk_1.item())
            risk_2 = w_2.T @ cov_ @ w_2
            risk_2 = np.sqrt(risk_2.item())
        elif rm == "MAD":
            risk_1 = MAD(a_1)
            risk_2 = MAD(a_2)
        elif rm == "GMD":
            risk_1 = GMD(a_1)
            risk_2 = GMD(a_2)
        elif rm == "MSV":
            risk_1 = SemiDeviation(a_1)
            risk_2 = SemiDeviation(a_2)
        elif rm == "FLPM":
            risk_1 = LPM(a_1, MAR=rf, p=1)
            risk_2 = LPM(a_2, MAR=rf, p=1)
        elif rm == "SLPM":
            risk_1 = LPM(a_1, MAR=rf, p=2)
            risk_2 = LPM(a_2, MAR=rf, p=2)
        elif rm == "VaR":
            risk_1 = VaR_Hist(a_1, alpha=alpha)
            risk_2 = VaR_Hist(a_2, alpha=alpha)
        elif rm == "CVaR":
            risk_1 = CVaR_Hist(a_1, alpha=alpha)
            risk_2 = CVaR_Hist(a_2, alpha=alpha)
        elif rm == "TG":
            risk_1 = TG(a_1, alpha=alpha, a_sim=a_sim)
            risk_2 = TG(a_2, alpha=alpha, a_sim=a_sim)
        elif rm == "EVaR":
            risk_1 = EVaR_Hist(a_1, alpha=alpha, solver=solver)[0]
            risk_2 = EVaR_Hist(a_2, alpha=alpha, solver=solver)[0]
        elif rm == "RLVaR":
            risk_1 = RLVaR_Hist(a_1, alpha=alpha, kappa=kappa, solver=solver)
            risk_2 = RLVaR_Hist(a_2, alpha=alpha, kappa=kappa, solver=solver)
        elif rm == "WR":
            risk_1 = WR(a_1)
            risk_2 = WR(a_2)
        elif rm == "VRG":
            risk_1 = VRG(a_1, alpha=alpha, beta=beta)
            risk_2 = VRG(a_2, alpha=alpha, beta=beta)
        elif rm == "CVRG":
            risk_1 = CVRG(a_1, alpha=alpha, beta=beta)
            risk_2 = CVRG(a_2, alpha=alpha, beta=beta)
        elif rm == "TGRG":
            risk_1 = TGRG(a_1, alpha=alpha, a_sim=a_sim, beta=beta, b_sim=b_sim)
            risk_2 = TGRG(a_2, alpha=alpha, a_sim=a_sim, beta=beta, b_sim=b_sim)
        elif rm == "EVRG":
            risk_1 = EVRG(a_1, alpha=alpha, beta=beta, solver=solver)
            risk_2 = EVRG(a_2, alpha=alpha, beta=beta, solver=solver)
        elif rm == "RVRG":
            risk_1 = RVRG(
                a_1, alpha=alpha, beta=beta, kappa=kappa, kappa_g=kappa_g, solver=solver
            )
            risk_2 = RVRG(
                a_2, alpha=alpha, beta=beta, kappa=kappa, kappa_g=kappa_g, solver=solver
            )
        elif rm == "RG":
            risk_1 = RG(a_1)
            risk_2 = RG(a_2)
        elif rm == "MDD":
            risk_1 = MDD_Abs(a_1)
            risk_2 = MDD_Abs(a_2)
        elif rm == "ADD":
            risk_1 = ADD_Abs(a_1)
            risk_2 = ADD_Abs(a_2)
        elif rm == "DaR":
            risk_1 = DaR_Abs(a_1, alpha=alpha)
            risk_2 = DaR_Abs(a_2, alpha=alpha)
        elif rm == "CDaR":
            risk_1 = CDaR_Abs(a_1, alpha=alpha)
            risk_2 = CDaR_Abs(a_2, alpha=alpha)
        elif rm == "EDaR":
            risk_1 = EDaR_Abs(a_1, alpha=alpha)[0]
            risk_2 = EDaR_Abs(a_2, alpha=alpha)[0]
        elif rm == "RLDaR":
            risk_1 = RLDaR_Abs(a_1, alpha=alpha, kappa=kappa, solver=solver)
            risk_2 = RLDaR_Abs(a_2, alpha=alpha, kappa=kappa, solver=solver)
        elif rm == "UCI":
            risk_1 = UCI_Abs(a_1)
            risk_2 = UCI_Abs(a_2)
        elif rm == "MDD_Rel":
            risk_1 = MDD_Rel(a_1)
            risk_2 = MDD_Rel(a_2)
        elif rm == "ADD_Rel":
            risk_1 = ADD_Rel(a_1)
            risk_2 = ADD_Rel(a_2)
        elif rm == "DaR_Rel":
            risk_1 = DaR_Rel(a_1, alpha=alpha)
            risk_2 = DaR_Rel(a_2, alpha=alpha)
        elif rm == "CDaR_Rel":
            risk_1 = CDaR_Rel(a_1, alpha=alpha)
            risk_2 = CDaR_Rel(a_2, alpha=alpha)
        elif rm == "EDaR_Rel":
            risk_1 = EDaR_Rel(a_1, alpha=alpha)[0]
            risk_2 = EDaR_Rel(a_2, alpha=alpha)[0]
        elif rm == "RLDaR_Rel":
            risk_1 = RLDaR_Rel(a_1, alpha=alpha, kappa=kappa, solver=solver)
            risk_2 = RLDaR_Rel(a_2, alpha=alpha, kappa=kappa, solver=solver)
        elif rm == "UCI_Rel":
            risk_1 = UCI_Rel(a_1)
            risk_2 = UCI_Rel(a_2)
        elif rm == "KT":
            risk_1 = Kurtosis(a_1) * 0.5
            risk_2 = Kurtosis(a_2) * 0.5
        elif rm == "SKT":
            risk_1 = SemiKurtosis(a_1) * 0.5
            risk_2 = SemiKurtosis(a_2) * 0.5

        RC_i = (risk_1 - risk_2) / (2 * d_i) * w_[i, 0]
        RC.append(RC_i)

    RC = np.array(RC, ndmin=1)

    return RC


def Risk_Margin(
    w,
    returns,
    cov=None,
    rm="MV",
    rf=0,
    alpha=0.05,
    a_sim=100,
    beta=None,
    b_sim=None,
    kappa=0.3,
    kappa_g=None,
    solver="CLARABEL",
):
    r"""
    Calculate the risk margin for each asset based on the risk measure
    selected.

    Parameters
    ----------
    w : DataFrame or Series of shape (n_assets, 1)
        Portfolio weights, where n_assets is the number of assets.
    returns : DataFrame or nd-array of shape (n_samples, n_features)
        Features matrix, where n_samples is the number of samples and
        n_features is the number of features.
    cov : DataFrame of shape (n_assets, n_assets)
        Covariance matrix, where n_assets is the number of assets.
    rm : str, optional
        Risk measure used in the denominator of the ratio. The default is
        'MV'. Possible values are:

        - 'MV': Standard Deviation.
        - 'KT': Square Root Kurtosis.
        - 'MAD': Mean Absolute Deviation.
        - 'GMD': Gini Mean Difference.
        - 'MSV': Semi Standard Deviation.
        - 'SKT': Square Root Semi Kurtosis.
        - 'FLPM': First Lower Partial Moment (Omega Ratio).
        - 'SLPM': Second Lower Partial Moment (Sortino Ratio).
        - 'VaR': Value at Risk.
        - 'CVaR': Conditional Value at Risk.
        - 'TG': Tail Gini.
        - 'EVaR': Entropic Value at Risk.
        - 'RLVaR': Relativistic Value at Risk. I recommend only use this function with MOSEK solver.
        - 'WR': Worst Realization (Minimax).
        - 'RG': Range of returns.
        - 'VRG' VaR range of returns.
        - 'CVRG': CVaR range of returns.
        - 'TGRG': Tail Gini range of returns.
        - 'EVRG': EVaR range of returns.
        - 'RVRG': RLVaR range of returns.
        - 'MDD': Maximum Drawdown of uncompounded cumulative returns (Calmar Ratio).
        - 'ADD': Average Drawdown of uncompounded cumulative returns.
        - 'DaR': Drawdown at Risk of uncompounded cumulative returns.
        - 'CDaR': Conditional Drawdown at Risk of uncompounded cumulative returns.
        - 'EDaR': Entropic Drawdown at Risk of uncompounded cumulative returns.
        - 'RLDaR': Relativistic Drawdown at Risk of uncompounded cumulative returns. I recommend only use this function with MOSEK solver.
        - 'UCI': Ulcer Index of uncompounded cumulative returns.
        - 'MDD_Rel': Maximum Drawdown of compounded cumulative returns (Calmar Ratio).
        - 'ADD_Rel': Average Drawdown of compounded cumulative returns.
        - 'CDaR_Rel': Conditional Drawdown at Risk of compounded cumulative returns.
        - 'EDaR_Rel': Entropic Drawdown at Risk of compounded cumulative returns.
        - 'RLDaR_Rel': Relativistic Drawdown at Risk of compounded cumulative returns. I recommend only use this function with MOSEK solver.
        - 'UCI_Rel': Ulcer Index of compounded cumulative returns.

    rf : float, optional
        Risk free rate. The default is 0.
    alpha : float, optional
        Significance level of VaR, CVaR, EVaR, RLVaR, DaR, CDaR, EDaR, RLDaR
        and Tail Gini of losses. The default is 0.05.
    a_sim : float, optional
        Number of CVaRs used to approximate Tail Gini of losses. The default is 100.
    beta : float, optional
        Significance level of VaR, CVaR, Tail Gini, EVaR and RLVaR of gains. If
        None it duplicates alpha value. The default is None.
    b_sim : float, optional
        Number of CVaRs used to approximate Tail Gini of gains. If None it
        duplicates a_sim value.
        The default is None.
    kappa : float, optional
        Deformation parameter of RLVaR and RLDaR for losses, must be between 0 and 1.
        The default is 0.3.
    kappa_g : float, optional
        Deformation parameter of RLVaR and RLDaR for gains, must be between 0 and 1.
        The default is None.
    solver: str, optional
        Solver available for CVXPY that supports exponential and power cone
        programming. Used to calculate EVaR, EVRG, EDaR, RLVaR, RVRG and RLDaR.
        The default value is None.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    Returns
    -------
    value : float
        Risk margin of the portfolio.

    """

    w_ = np.array(w, ndmin=2)
    if w_.shape[0] == 1 and w_.shape[1] > 1:
        w_ = w_.T
    if w_.shape[0] > 1 and w_.shape[1] > 1:
        raise ValueError("weights must have n_assets x 1 size")

    if isinstance(returns, pd.Series):
        returns_ = returns.to_frame()
    elif isinstance(returns, pd.DataFrame):
        returns_ = returns.to_numpy()
    else:
        returns_ = np.array(returns, ndmin=2)

    if cov is None:
        cov_ = np.array(np.cov(returns_, rowvar=False), ndmin=2)
    else:
        cov_ = np.array(cov, ndmin=2)

    RM = []
    if rm in ["RLVaR", "RLDaR"]:
        d_i = 0.0001
    else:
        d_i = 0.0000001

    for i in range(0, w_.shape[0]):
        delta = np.zeros((w_.shape[0], 1))
        delta[i, 0] = d_i
        w_1 = w_ + delta
        w_2 = w_ - delta
        a_1 = returns_ @ w_1
        a_2 = returns_ @ w_2
        if rm == "MV":
            risk_1 = w_1.T @ cov_ @ w_1
            risk_1 = np.sqrt(risk_1.item())
            risk_2 = w_2.T @ cov_ @ w_2
            risk_2 = np.sqrt(risk_2.item())
        elif rm == "MAD":
            risk_1 = MAD(a_1)
            risk_2 = MAD(a_2)
        elif rm == "GMD":
            risk_1 = GMD(a_1)
            risk_2 = GMD(a_2)
        elif rm == "MSV":
            risk_1 = SemiDeviation(a_1)
            risk_2 = SemiDeviation(a_2)
        elif rm == "FLPM":
            risk_1 = LPM(a_1, MAR=rf, p=1)
            risk_2 = LPM(a_2, MAR=rf, p=1)
        elif rm == "SLPM":
            risk_1 = LPM(a_1, MAR=rf, p=2)
            risk_2 = LPM(a_2, MAR=rf, p=2)
        elif rm == "VaR":
            risk_1 = VaR_Hist(a_1, alpha=alpha)
            risk_2 = VaR_Hist(a_2, alpha=alpha)
        elif rm == "CVaR":
            risk_1 = CVaR_Hist(a_1, alpha=alpha)
            risk_2 = CVaR_Hist(a_2, alpha=alpha)
        elif rm == "TG":
            risk_1 = TG(a_1, alpha=alpha, a_sim=a_sim)
            risk_2 = TG(a_2, alpha=alpha, a_sim=a_sim)
        elif rm == "EVaR":
            risk_1 = EVaR_Hist(a_1, alpha=alpha, solver=solver)[0]
            risk_2 = EVaR_Hist(a_2, alpha=alpha, solver=solver)[0]
        elif rm == "RLVaR":
            risk_1 = RLVaR_Hist(a_1, alpha=alpha, kappa=kappa, solver=solver)
            risk_2 = RLVaR_Hist(a_2, alpha=alpha, kappa=kappa, solver=solver)
        elif rm == "WR":
            risk_1 = WR(a_1)
            risk_2 = WR(a_2)
        elif rm == "VRG":
            risk_1 = VRG(a_1, alpha=alpha, beta=beta)
            risk_2 = VRG(a_2, alpha=alpha, beta=beta)
        elif rm == "CVRG":
            risk_1 = CVRG(a_1, alpha=alpha, beta=beta)
            risk_2 = CVRG(a_2, alpha=alpha, beta=beta)
        elif rm == "TGRG":
            risk_1 = TGRG(a_1, alpha=alpha, a_sim=a_sim, beta=beta, b_sim=b_sim)
            risk_2 = TGRG(a_2, alpha=alpha, a_sim=a_sim, beta=beta, b_sim=b_sim)
        elif rm == "EVRG":
            risk_1 = EVRG(a_1, alpha=alpha, beta=beta, solver=solver)
            risk_2 = EVRG(a_2, alpha=alpha, beta=beta, solver=solver)
        elif rm == "RVRG":
            risk_1 = RVRG(
                a_1, alpha=alpha, beta=beta, kappa=kappa, kappa_g=kappa_g, solver=solver
            )
            risk_2 = RVRG(
                a_2, alpha=alpha, beta=beta, kappa=kappa, kappa_g=kappa_g, solver=solver
            )
        elif rm == "RG":
            risk_1 = RG(a_1)
            risk_2 = RG(a_2)
        elif rm == "MDD":
            risk_1 = MDD_Abs(a_1)
            risk_2 = MDD_Abs(a_2)
        elif rm == "ADD":
            risk_1 = ADD_Abs(a_1)
            risk_2 = ADD_Abs(a_2)
        elif rm == "DaR":
            risk_1 = DaR_Abs(a_1, alpha=alpha)
            risk_2 = DaR_Abs(a_2, alpha=alpha)
        elif rm == "CDaR":
            risk_1 = CDaR_Abs(a_1, alpha=alpha)
            risk_2 = CDaR_Abs(a_2, alpha=alpha)
        elif rm == "EDaR":
            risk_1 = EDaR_Abs(a_1, alpha=alpha, solver=solver)[0]
            risk_2 = EDaR_Abs(a_2, alpha=alpha, solver=solver)[0]
        elif rm == "RLDaR":
            risk_1 = RLDaR_Abs(a_1, alpha=alpha, kappa=kappa, solver=solver)
            risk_2 = RLDaR_Abs(a_2, alpha=alpha, kappa=kappa, solver=solver)
        elif rm == "UCI":
            risk_1 = UCI_Abs(a_1)
            risk_2 = UCI_Abs(a_2)
        elif rm == "MDD_Rel":
            risk_1 = MDD_Rel(a_1)
            risk_2 = MDD_Rel(a_2)
        elif rm == "ADD_Rel":
            risk_1 = ADD_Rel(a_1)
            risk_2 = ADD_Rel(a_2)
        elif rm == "DaR_Rel":
            risk_1 = DaR_Rel(a_1, alpha=alpha)
            risk_2 = DaR_Rel(a_2, alpha=alpha)
        elif rm == "CDaR_Rel":
            risk_1 = CDaR_Rel(a_1, alpha=alpha)
            risk_2 = CDaR_Rel(a_2, alpha=alpha)
        elif rm == "EDaR_Rel":
            risk_1 = EDaR_Rel(a_1, alpha=alpha, solver=solver)[0]
            risk_2 = EDaR_Rel(a_2, alpha=alpha, solver=solver)[0]
        elif rm == "RLDaR_Rel":
            risk_1 = RLDaR_Rel(a_1, alpha=alpha, kappa=kappa, solver=solver)
            risk_2 = RLDaR_Rel(a_2, alpha=alpha, kappa=kappa, solver=solver)
        elif rm == "UCI_Rel":
            risk_1 = UCI_Rel(a_1)
            risk_2 = UCI_Rel(a_2)
        elif rm == "KT":
            risk_1 = Kurtosis(a_1) * 0.5
            risk_2 = Kurtosis(a_2) * 0.5
        elif rm == "SKT":
            risk_1 = SemiKurtosis(a_1) * 0.5
            risk_2 = SemiKurtosis(a_2) * 0.5

        RM_i = (risk_1 - risk_2) / (2 * d_i)
        RM.append(RM_i)

    RM = np.array(RM, ndmin=1)

    return RM


def Factors_Risk_Contribution(
    w,
    returns,
    factors,
    cov=None,
    B=None,
    const=False,
    rm="MV",
    rf=0,
    alpha=0.05,
    a_sim=100,
    beta=None,
    b_sim=None,
    kappa=0.3,
    kappa_g=None,
    solver="CLARABEL",
    feature_selection="stepwise",
    stepwise="Forward",
    criterion="pvalue",
    threshold=0.05,
    n_components=0.95,
):
    r"""
    Calculate the risk contribution for each factor based on the selected risk measure.

    Parameters
    ----------
    w : DataFrame or Series of shape (n_assets, 1)
        Portfolio weights, where n_assets is the number of assets.
    returns : DataFrame or nd-array of shape (n_samples, n_features)
        Features matrix, where n_samples is the number of samples and
        n_features is the number of features.
    factors : DataFrame or nd-array of shape (n_samples, n_factors)
        Factors matrix, where n_samples is the number of samples and
        n_factors is the number of factors.
    cov : DataFrame of shape (n_assets, n_assets)
        Covariance matrix, where n_assets is the number of assets.
    B : DataFrame of shape (n_assets, n_factors), optional
        Loadings matrix, where n_assets is the number assets and n_factors is
        the number of risk factors. If is not specified, is estimated using
        stepwise regression. The default is None.
    const : bool, optional
        Indicate if the loadings matrix has a constant.
        The default is False.
    rm : str, optional
        Risk measure used in the denominator of the ratio. The default is
        'MV'. Possible values are:

        - 'MV': Standard Deviation.
        - 'KT': Square Root Kurtosis.
        - 'MAD': Mean Absolute Deviation.
        - 'GMD': Gini Mean Difference.
        - 'MSV': Semi Standard Deviation.
        - 'SKT': Square Root Semi Kurtosis.
        - 'FLPM': First Lower Partial Moment (Omega Ratio).
        - 'SLPM': Second Lower Partial Moment (Sortino Ratio).
        - 'VaR': Value at Risk.
        - 'CVaR': Conditional Value at Risk.
        - 'TG': Tail Gini.
        - 'EVaR': Entropic Value at Risk.
        - 'RLVaR': Relativistic Value at Risk. I recommend only use this function with MOSEK solver.
        - 'WR': Worst Realization (Minimax).
        - 'RG': Range of returns.
        - 'VRG' VaR range of returns.
        - 'CVRG': CVaR range of returns.
        - 'TGRG': Tail Gini range of returns.
        - 'EVRG': EVaR range of returns.
        - 'RVRG': RLVaR range of returns. I recommend only use this function with MOSEK solver.
        - 'MDD': Maximum Drawdown of uncompounded cumulative returns (Calmar Ratio).
        - 'ADD': Average Drawdown of uncompounded cumulative returns.
        - 'DaR': Drawdown at Risk of uncompounded cumulative returns.
        - 'CDaR': Conditional Drawdown at Risk of uncompounded cumulative returns.
        - 'EDaR': Entropic Drawdown at Risk of uncompounded cumulative returns.
        - 'RLDaR': Relativistic Drawdown at Risk of uncompounded cumulative returns. I recommend only use this function with MOSEK solver.
        - 'UCI': Ulcer Index of uncompounded cumulative returns.
        - 'MDD_Rel': Maximum Drawdown of compounded cumulative returns (Calmar Ratio).
        - 'ADD_Rel': Average Drawdown of compounded cumulative returns.
        - 'CDaR_Rel': Conditional Drawdown at Risk of compounded cumulative returns.
        - 'EDaR_Rel': Entropic Drawdown at Risk of compounded cumulative returns.
        - 'RLDaR_Rel': Relativistic Drawdown at Risk of compounded cumulative returns. I recommend only use this function with MOSEK solver.
        - 'UCI_Rel': Ulcer Index of compounded cumulative returns.

    rf : float, optional
        Risk free rate. The default is 0.
    alpha : float, optional
        Significance level of VaR, CVaR, EVaR, RLVaR, DaR, CDaR, EDaR, RLDaR and Tail Gini of losses.
        The default is 0.05.
    a_sim : float, optional
        Number of CVaRs used to approximate Tail Gini of losses. The default is 100.
    beta : float, optional
        Significance level of CVaR and Tail Gini of gains. If None it duplicates alpha value.
        The default is None.
    b_sim : float, optional
        Number of CVaRs used to approximate Tail Gini of gains. If None it duplicates a_sim value.
        The default is None.
    kappa : float, optional
        Deformation parameter of RLVaR and RLDaR for losses, must be between 0 and 1.
        The default is 0.3.
    kappa_g : float, optional
        Deformation parameter of RLVaR and RLDaR for gains, must be between 0 and 1.
        The default is None.
    solver: str, optional
        Solver available for CVXPY that supports power cone programming. Used to calculate RLVaR and RLDaR.
        The default value is None.
    feature_selection: str 'stepwise' or 'PCR', optional
        Indicate the method used to estimate the loadings matrix.
        The default is 'stepwise'.
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

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    Returns
    -------
    value : float
        Risk measure of the portfolio.

    """
    w_ = np.array(w, ndmin=2)
    if w_.shape[0] == 1 and w_.shape[1] > 1:
        w_ = w_.T
    if w_.shape[0] > 1 and w_.shape[1] > 1:
        raise ValueError("weights must have n_assets x 1 size")

    if returns.index.tolist() != factors.index.tolist():
        raise ValueError("returns and factors must have same dates.")

    RM = Risk_Margin(
        w=w_,
        returns=returns,
        cov=cov,
        rm=rm,
        rf=rf,
        alpha=alpha,
        a_sim=a_sim,
        beta=beta,
        b_sim=b_sim,
        kappa=kappa,
        solver=solver,
    ).reshape(-1, 1)

    if B is None:
        B = pe.loadings_matrix(
            X=factors,
            Y=returns,
            feature_selection=feature_selection,
            stepwise=stepwise,
            criterion=criterion,
            threshold=threshold,
            n_components=n_components,
        )
        const = True
    elif not isinstance(B, pd.DataFrame):
        raise ValueError("B must be a DataFrame")

    if const == True or factors.shape[1] + 1 == B.shape[1]:
        B = B.iloc[:, 1:].to_numpy()

    if feature_selection == "PCR":
        scaler = StandardScaler()
        scaler.fit(factors)
        factors_std = scaler.transform(factors)
        if n_components > 0 and n_components < 1:
            pca = PCA(n_components=n_components)
        elif n_components >= 1:
            pca = PCA(n_components=int(n_components))
        pca.fit(factors_std)
        V_p = pca.components_.T
        std = np.array(np.std(factors, axis=0, ddof=1), ndmin=2)
        B = (pinv(V_p) @ (B.T * std.T)).T

    B1 = pinv(B.T)
    B2 = pinv(null_space(B.T).T)
    B3 = pinv(B2.T)

    RC_F = (B.T @ w_) * (B1.T @ RM)
    RC_OF = np.array(((B2.T @ w_) * (B3.T @ RM)).sum(), ndmin=2)
    RC_F = np.vstack([RC_F, RC_OF]).ravel()

    return RC_F


def BrinsonAttribution(
    prices,
    w,
    wb,
    start,
    end,
    asset_classes,
    classes_col,
    method="nearest",
):
    r"""
    Creates a DataFrame with the Brinson Performance Attribution per class and
    aggregate based on :cite:`f-Brinson1985`.

    Parameters
    ----------
    prices : DataFrame of shape (n_samples, n_assets)
        Assets prices DataFrame, where n_samples is the number of
        observations and n_assets is the number of assets.
    w : DataFrame  or Series of shape (n_assets, 1)
        A portfolio specified by the user.
    wb : DataFrame  or Series of shape (n_assets, 1)
        A benchmark specified by the user.
    start : str
        Start date in format 'YYYY-MM-DD' specified by the user.
    end : str
        End date in format 'YYYY-MM-DD' specified by the user.
    asset_classes : DataFrame of shape (n_assets, n_cols)
        Asset's classes DataFrame, where n_assets is the number of assets and
        n_cols is the number of columns of the DataFrame where the first column
        is the asset list and the next columns are the different asset's
        classes sets. It is only used when kind value is 'classes'. The default
        value is None.
    classes_col : str or int
        If value is str, it is the column name of the set of classes from
        asset_classes dataframe. If value is int, it is the column number of
        the set of classes from asset_classes dataframe. The default
        value is None.
    method : str
        Method used to calculate the nearest start or end dates in case one of
        them is not in prices DataFrame. The default value is 'nearest'.
        See `get_indexer <https://pandas.pydata.org/docs/reference/api/pandas.Index.get_indexer.html#pandas.Index.get_indexer>`__ for more details.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    Returns
    -------
    BrinAttr : DataFrame
        A DataFrame with the Brinson Performance Attribution per class and aggregate.

    (start_, end_) : tuple
        Start and end dates calculated using get_indexer method in string format.


    Example
    -------
    ::

        BrinAttr, (start, end) = BrinsonAttribution(
            prices=data,
            w=w,
            wb=wb,
            start='2019-01-07',
            end='2019-12-06',
            asset_classes=asset_classes,
            classes_col='Industry',
            )

    .. image:: images/BrinAttr.png


    """

    if not isinstance(prices, pd.DataFrame):
        raise ValueError("prices must be a DataFrame")

    if not isinstance(w, pd.DataFrame):
        if isinstance(w, pd.Series):
            wp_ = w.to_frame()
        else:
            raise ValueError("w must be a one column DataFrame or Series")
    else:
        if w.shape[0] == 1:
            wp_ = w.T.copy()
        elif w.shape[1] == 1:
            wp_ = w.copy()
        else:
            raise ValueError("w must be a one column DataFrame or Series")

    if not isinstance(wb, pd.DataFrame):
        if isinstance(wb, pd.Series):
            wb_ = wb.to_frame()
        else:
            raise ValueError("w must be a one column DataFrame or Series")
    else:
        if wb.shape[0] == 1:
            wb_ = wb.T.copy()
        elif wb.shape[1] == 1:
            wb_ = wb.copy()
        else:
            raise ValueError("w must be a one column DataFrame or Series")

    if not isinstance(asset_classes, pd.DataFrame):
        raise ValueError("asset_classes must be a DataFrame")
    else:
        if asset_classes.shape[1] < 2:
            raise ValueError("asset_classes must have at least two columns")
        classes = asset_classes.columns.tolist()
        if isinstance(classes_col, str) and classes_col in classes:
            col = classes_col
        elif isinstance(classes_col, int) and classes[classes_col] in classes:
            col = classes[classes_col]
        else:
            raise ValueError(
                "classes_col must be a valid column or column position of asset_classes"
            )

    prices_ = prices.copy()
    prices_.index = prices_.index.tz_localize(None)

    start_ = prices_.index.get_indexer([pd.Timestamp(start)], method=method)
    end_ = prices_.index.get_indexer([pd.Timestamp(end)], method=method)

    p1 = prices_.iloc[start_].to_numpy().reshape(-1, 1)
    p2 = prices_.iloc[end_].to_numpy().reshape(-1, 1)
    p3 = p2 / p1 - 1

    wp_ = wp_.to_numpy().reshape(-1, 1)
    wb_ = wb_.to_numpy().reshape(-1, 1)

    Rp = (p3.T @ wp_).item()
    Rb = (p3.T @ wb_).item()

    classes = asset_classes[col].tolist()
    unique_classes = list(set(classes))
    unique_classes.sort()

    labels = [
        "Asset Allocation",
        "Security Selection",
        "Interaction",
        "Total Excess Return",
    ]
    BrinAttr = pd.DataFrame([], index=labels)

    for i in unique_classes:
        sets_i = []
        for j in classes:
            sets_i.append(i == j)
        sets_i = np.array(sets_i, dtype=int).reshape(-1, 1)

        wb_i = (sets_i.T @ wb_).item()
        wp_i = (sets_i.T @ wp_).item()

        Rb_i = (np.multiply(p3, sets_i).T @ wb_).item() / wb_i
        Rp_i = (np.multiply(p3, sets_i).T @ wp_).item() / wp_i

        AAE_i = (wp_i - wb_i) * (Rb_i - Rb)
        SSE_i = wb_i * (Rp_i - Rb_i)
        IE_i = (wp_i - wb_i) * (Rp_i - Rb_i)
        TER_i = AAE_i + SSE_i + IE_i

        BrinAttr_i = pd.DataFrame(
            [AAE_i, SSE_i, IE_i, TER_i], index=labels, columns=[i]
        )

        BrinAttr = pd.concat([BrinAttr, BrinAttr_i], axis=1)

    total = BrinAttr.sum(axis=1).to_frame()
    total.columns = ["Total"]

    BrinAttr = pd.concat([BrinAttr, total], axis=1)

    start_ = prices_.index.tolist()[start_.item()].strftime("%Y-%m-%d")
    end_ = prices_.index.tolist()[end_.item()].strftime("%Y-%m-%d")

    return BrinAttr, (start_, end_)
