import numpy as np
from scipy.optimize import minimize
from scipy.optimize import Bounds


__all__ = [
    "MAD",
    "SemiDeviation",
    "VaR_Hist",
    "CVaR_Hist",
    "WR",
    "LPM",
    "Entropic_RM",
    "EVaR_Hist",
    "MaxAbsDD",
    "AvgAbsDD",
    "ConAbsDD",
    "MaxRelDD",
    "AvgRelDD",
    "ConRelDD",
    "Sharpe_Risk",
    "Sharpe",
    "Risk_Contribution",
]


def MAD(X):
    r"""
    Calculates the Mean Absolute Deviation (MAD) of a returns series.

    .. math::
        \text{MAD}(X) = \frac{1}{T}\sum_{t=1}^{T}
        | X_{t} - \mathbb{E}(X_{t}) |

    Parameters
    ----------
    X : 1d-array 
        a returns series, must have Tx1 size.

    Returns
    -------
    value : float    
        MAD of a returns series.
        
    Raises
    ------
    ValueError
        When the value cannot be calculated.
        
    Examples
    --------
    Examples should be written in doctest format, and should illustrate how
    to use the function.

    >>> print([i for i in example_generator(4)])
    [0, 1, 2, 3]
      
    """

    a = np.array(X, ndmin=2)
    if a.shape[0] == 1 and a.shape[1] > 1:
        a = a.T
    if a.shape[0] > 1 and a.shape[1] > 1:
        raise ValueError("returns must have Tx1 size")

    value = np.mean(np.absolute(a - np.mean(a, axis=0)), axis=0)
    value = value.item()

    return value


def SemiDeviation(X):
    r"""
    Calculates the Semi Deviation of a returns series.

    .. math::
        \text{SemiDev}(X) = \left [ \frac{1}{T-1}\sum_{t=1}^{T}
        (X_{t} - \mathbb{E}(X_{t}))^2 \right ]^{1/2}    

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

    mu = np.mean(a, axis=0)
    value = mu - a
    n = value.shape[0] - 1
    value = np.sum(np.power(value[np.where(value <= mu)], 2)) / n
    value = np.power(value, 0.5).item()

    return value


def VaR_Hist(X, alpha=0.01):
    r"""
    Calculates the Value at Risk (VaR) of a returns series.

    .. math::
        \text{VaR}_{\alpha}(X) = -\inf_{t \in (0,T)} \left \{ X_{t} \in
        \mathbb{R}: F_{X}(X_{t})>\alpha \right \}

    Parameters
    ----------
    X : 1d-array 
        Returns series, must have Tx1 size.
    alpha : float, optional
        Significance level of VaR. The default is 0.01.
                        
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
    value = value.item()

    return value


def CVaR_Hist(X, alpha=0.01):
    r"""
    Calculates the Conditional Value at Risk (CVaR) of a returns series.

    .. math::
        \text{CVaR}_{\alpha}(X) = \text{VaR}_{\alpha}(X) +
        \frac{1}{\alpha T} \sum_{t=1}^{T} \max(-X_{t} -
        \text{VaR}_{\alpha}(X), 0)

    Parameters
    ----------
    X : 1d-array 
        Returns series, must have Tx1 size.
    alpha : float, optional
        Significance level of CVaR. The default is 0.01.
                        
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
    value = value.item()

    return value


def WR(X):
    r"""
    Calculates the Worst Realization (WR) or Worst Scenario of a returns series.

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
    value = value.item()

    return value


def LPM(X, MAR=0, p=1):
    r"""
    Calculates the p-th Lower Partial Moment of a returns series.

    .. math::
        \text{LPM}(X, \text{MAR}, p) = \left [ \frac{1}{T}\sum_{t=1}^{T}
        \max(\text{MAR} - X_{t}, 0) \right ]^{\frac{1}{p}}    

    Where:
    
    :math:`\text{MAR}` is the minimum acceptable return.

    Parameters
    ----------
    X : 1d-array 
        Returns series, must have Tx1 size.
    MAR : float, optional
        Minimum acceptable return. The default is 0.
    p : float, optional
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

    value = MAR - a

    if p > 1:
        n = value.shape[0] - 1
    else:
        n = value.shape[0]

    value = np.sum(np.power(value[np.where(value > 0)], p)) / n
    value = np.power(value, 1 / p).item()

    return value


def Entropic_RM(X, theta=1):
    r"""
    Calculates the Entropic Risk Measure (ERM) of a returns series.

    .. math::
        \text{ERM}(X) = \theta \log\left(\mathbb{E}
        [e^{-\frac{1}{\theta} X}]\right)        
        
    Parameters
    ----------
    X : 1d-array 
        Returns series, must have Tx1 size.
    theta : float, optional
        Risk aversion parameter, must be greater than zero. The default is 1.
                        
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

    value = np.mean(np.exp(-1 / theta * np.array(a)), axis=0)
    value = theta * (np.log(value))
    value = value.item()

    return value


def _Entropic_RM(X, theta=1, alpha=0.01):
    a = np.array(X, ndmin=2)
    if a.shape[0] == 1 and a.shape[1] > 1:
        a = a.T
    if a.shape[0] > 1 and a.shape[1] > 1:
        raise ValueError("returns must have Tx1 size")

    value = np.mean(np.exp(-1 / theta * np.array(a)), axis=0)
    value = theta * (np.log(value) - np.log(alpha))
    value = value.item()

    return value


def EVaR_Hist(X, alpha=0.01):
    r"""
    Calculates the Entropic Value at Risk (EVaR) of a returns series.

    .. math::
        \text{EVaR}_{\alpha}(X) = \inf_{z>0} \left \{ z^{-1}
        \ln \left (\frac{M_X(z)}{\alpha} \right ) \right \}
    
    Where:
    
    :math:`M_X(z)` is the moment generating function of X.
    
    Parameters
    ----------
    X : 1d-array 
        Returns series, must have Tx1 size.
    alpha : float, optional
        Significance level of EVaR. The default is 0.01.
                        
    Raises
    ------
    ValueError
        When the value cannot be calculated.

    Returns
    -------
    value : float
        EVaR of a returns series.
    
    """

    a = np.array(X, ndmin=2)
    if a.shape[0] == 1 and a.shape[1] > 1:
        a = a.T
    if a.shape[0] > 1 and a.shape[1] > 1:
        raise ValueError("returns must have Tx1 size")

    bnd = Bounds([0.00000000001], [np.inf])
    result = minimize(_Entropic_RM, [0.01], args=(X, alpha), bounds=bnd)
    t = result.x
    t = t.item()
    value = _Entropic_RM(t, X, alpha)
    return value


def MaxAbsDD(X):
    r"""
    Calculates the Maximum Drawdown (MDD) of a returns series
    using uncumpound cumulated returns.

    .. math::
        \text{MDD}(X) = \max_{j \in (0,T)} \left [\max_{t \in (0,T)}
        \left ( \sum_{i=0}^{t}X_{i} - \sum_{i=0}^{j}X_{i} \right ) \right ]
    
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
        MDD of a uncumpound cumulated returns.
    
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

    value = value.item()

    return value


def AvgAbsDD(X):
    r"""
    Calculates the Average Drawdown (ADD) of a returns series
    using uncumpound cumulated returns.

    .. math::
        \text{ADD}(X) = \frac{1}{T}\sum_{i=0}^{T}\max_{t \in (0,T)}
        \left ( \sum_{i=0}^{t}X_{i} - \sum_{i=0}^{j}X_{i} \right ) 
    
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
        ADD of a uncumpound cumulated returns.
    
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
        value = value / n

    value = value.item()

    return value


def ConAbsDD(X, alpha=0.01):
    r"""
    Calculates the Conditional Drawdown at Risk (CDaR) of a returns series
    using uncumpound cumulated returns.

    .. math::
        \text{CDaR}_{\alpha}(X) = \text{DaR}_{\alpha}(X) + \frac{1}{\alpha T}
        \sum_{i=0}^{T} \max \left [ \max_{t \in (0,T)}
        \left ( \sum_{i=0}^{t}X_{i} - \sum_{i=0}^{j}X_{i} \right )
        - \text{DaR}_{\alpha}(X), 0 \right ]
    
    Where:

    :math:`\text{DaR}_{\alpha}` is the Drawdown at Risk of an uncumpound
    cumulated return series :math:`X`.  

    Parameters
    ----------
    X : 1d-array
        Returns series, must have Tx1 size..
    alpha : float, optional
        Significance level of CDaR. The default is 0.01.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    Returns
    -------
    value : float
        CDaR of a uncumpound cumulated returns series.

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
    value = value.item()

    return value


def MaxRelDD(X):
    r"""
    Calculates the Maximum Drawdown (MDD) of a returns series
    using cumpound cumulated returns.

    .. math::
        \text{MDD}(X) = \max_{j \in (0,T)}\left[\max_{t \in (0,T)}
        \left ( \prod_{i=0}^{t}(1+X_{i}) - \prod_{i=0}^{j}(1+X_{i}) \right ) \right]
    
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
        MDD of a cumpound cumulated returns.

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

    value = value.item()

    return value


def AvgRelDD(X):
    r"""
    Calculates the Average Drawdown (ADD) of a returns series
    using cumpound acumulated returns.

    .. math::
        \text{ADD}(X) = \frac{1}{T}\sum_{i=0}^{T}\max_{t \in (0,T)}
        \left ( \prod_{i=0}^{t}(1+X_{i}) - \prod_{i=0}^{j}(1+X_{i}) \right ) 
    
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
        ADD of a cumpound acumulated returns.
    
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
        value = value / n

    value = value.item()

    return value


def ConRelDD(X, alpha=0.01):
    r"""
    Calculates the Conditional Drawdown at Risk (CDaR) of a returns series
    using cumpound cumulated returns.

    .. math::
        \text{CDaR}_{\alpha}(X) = \text{DaR}_{\alpha}(X) + \frac{1}{\alpha T}
        \sum_{i=0}^{T} \max \left [ \max_{t \in (0,T)}
        \left ( \prod_{i=0}^{t}(1+X_{i}) - \prod_{i=0}^{j}(1+X_{i}) \right )
        - \text{DaR}_{\alpha}(X), 0 \right ]
    
    Where:

    :math:`\text{DaR}_{\alpha}` is the Drawdown at Risk of a cumpound
    acumulated return series :math:`X`.       

    Parameters
    ----------
    X : 1d-array
        Returns series, must have Tx1 size..
    alpha : float, optional
        Significance level of CDaR. The default is 0.01.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    Returns
    -------
    value : float
        CDaR of a cumpound cumulated returns series.

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
    value = value.item()

    return value


###############################################################################
# Risk Adjusted Return Ratios
###############################################################################


def Sharpe_Risk(w, cov=None, returns=None, rm="MV", rf=0, alpha=0.01):
    r"""
    Calculate the risk measure available on the Sharpe function.

    Parameters
    ----------
    w : DataFrame or 1d-array of shape (n_assets, 1)
        Weights matrix, where n_assets is the number of assets.
    cov : DataFrame or nd-array of shape (n_features, n_features)
        Covariance matrix, where n_features is the number of features.
    returns : DataFrame or nd-array of shape (n_samples, n_features)
        Features matrix, where n_samples is the number of samples and 
        n_features is the number of features.
    rm : str, optional
        Risk measure used in the denominator of the ratio. The default is
        'MV'. Posible values are:
            
        - 'MV': Standard Deviation.
        - 'MAD': Mean Absolute Deviation.
        - 'MSV': Semi Standard Deviation.
        - 'FLPM': First Lower Partial Moment (Omega Ratio).
        - 'SLPM': Second Lower Partial Moment (Sortino Ratio).
        - 'VaR': Value at Risk.
        - 'CVaR': Conditional Value at Risk.
        - 'WR': Worst Realization (Minimax)
        - 'MDD': Maximum Drawdown of uncompounded returns (Calmar Ratio).
        - 'ADD': Average Drawdown of uncompounded returns.
        - 'CDaR': Conditional Drawdown at Risk of uncompounded returns.
        
    rf : float, optional
        Risk free rate. The default is 0.
    **kwargs : dict
        Other arguments that depends on the risk measure.

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
    if cov is not None:
        cov_ = np.array(cov, ndmin=2)
    if returns is not None:
        returns_ = np.array(returns, ndmin=2)

    a = returns_ @ w_
    if rm == "MV":
        risk = w_.T @ cov_ @ w_
        risk = np.sqrt(risk.item())
    elif rm == "MAD":
        risk = MAD(a)
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
    elif rm == "WR":
        risk = WR(a)
    elif rm == "MDD":
        risk = MaxAbsDD(a)
    elif rm == "ADD":
        risk = AvgAbsDD(a)
    elif rm == "CDaR":
        risk = ConAbsDD(a, alpha=alpha)

    value = risk

    return value


def Sharpe(w, mu, cov=None, returns=None, rm="MV", rf=0, alpha=0.01):
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
    
    w : DataFrame or 1d-array of shape (n_assets, 1)
        Weights matrix, where n_assets is the number of assets.
    mu : DataFrame or nd-array of shape (1, n_assets)
        Vector of expected returns, where n_assets is the number of assets.
    cov : DataFrame or nd-array of shape (n_features, n_features)
        Covariance matrix, where n_features is the number of features.
    returns : DataFrame or nd-array of shape (n_samples, n_features)
        Features matrix, where n_samples is the number of samples and 
        n_features is the number of features.
    rm : str, optional
        Risk measure used in the denominator of the ratio. The default is
        'MV'. Posible values are:
            
        - 'MV': Standard Deviation.
        - 'MAD': Mean Absolute Deviation.
        - 'MSV': Semi Standard Deviation.
        - 'FLPM': First Lower Partial Moment (Omega Ratio).
        - 'SLPM': Second Lower Partial Moment (Sortino Ratio).
        - 'VaR': Value at Risk.
        - 'CVaR': Conditional Value at Risk.
        - 'WR': Worst Realization (Minimax)
        - 'MDD': Maximum Drawdown of uncompounded returns (Calmar Ratio).
        - 'ADD': Average Drawdown of uncompounded returns.
        - 'CDaR': Conditional Drawdown at Risk of uncompounded returns.
        
    rf : float, optional
        Risk free rate. The default is 0.
    **kwargs : dict
        Other arguments that depends on the risk measure.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    Returns
    -------
    value : float
        Risk adjusted return ratio of :math:`X`.

    """

    if cov is None and rm == "MV":
        raise ValueError("covariance matrix is necessary to calculate the sharpe ratio")
    elif returns is None and rm != "MV":
        raise ValueError(
            "returns scenarios are necessary to calculate the sharpe ratio"
        )

    w_ = np.array(w, ndmin=2)
    mu_ = np.array(mu, ndmin=2)

    if cov is not None:
        cov_ = np.array(cov, ndmin=2)
    if returns is not None:
        returns_ = np.array(returns, ndmin=2)

    ret = mu_ @ w_
    ret = ret.item()

    risk = Sharpe_Risk(w, cov=cov_, returns=returns_, rm=rm, rf=rf, alpha=alpha)

    value = (ret - rf) / risk

    return value


###############################################################################
# Risk Contribution Vectors
###############################################################################


def Risk_Contribution(w, cov=None, returns=None, rm="MV", rf=0, alpha=0.01):
    r"""
    Calculate the risk contribution for each asset based on the risk measure
    selected.

    Parameters
    ----------
    w : DataFrame or 1d-array of shape (n_assets, 1)
        Weights matrix, where n_assets is the number of assets.
    cov : DataFrame or nd-array of shape (n_features, n_features)
        Covariance matrix, where n_features is the number of features.
    returns : DataFrame or nd-array of shape (n_samples, n_features)
        Features matrix, where n_samples is the number of samples and 
        n_features is the number of features.
    rm : str, optional
        Risk measure used in the denominator of the ratio. The default is
        'MV'. Posible values are:
            
        - 'MV': Standard Deviation.
        - 'MAD': Mean Absolute Deviation.
        - 'MSV': Semi Standard Deviation.
        - 'FLPM': First Lower Partial Moment (Omega Ratio).
        - 'SLPM': Second Lower Partial Moment (Sortino Ratio).
        - 'VaR': Value at Risk.
        - 'CVaR': Conditional Value at Risk.
        - 'WR': Worst Realization (Minimax)
        - 'MDD': Maximum Drawdown of uncompounded returns (Calmar Ratio).
        - 'ADD': Average Drawdown of uncompounded returns.
        - 'CDaR': Conditional Drawdown at Risk of uncompounded returns.
        
    rf : float, optional
        Risk free rate. The default is 0.
    **kwargs : dict
        Other arguments that depends on the risk measure.

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
    if cov is not None:
        cov_ = np.array(cov, ndmin=2)
    if returns is not None:
        returns_ = np.array(returns, ndmin=2)

    # risk = Sharpe_Risk(w, cov=cov_, returns=returns_, rm=rm, rf=rf, alpha=alpha)

    RC = []
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
        elif rm == "WR":
            risk_1 = WR(a_1)
            risk_2 = WR(a_2)
        elif rm == "MDD":
            risk_1 = MaxAbsDD(a_1)
            risk_2 = MaxAbsDD(a_2)
        elif rm == "ADD":
            risk_1 = AvgAbsDD(a_1)
            risk_2 = AvgAbsDD(a_2)
        elif rm == "CDaR":
            risk_1 = ConAbsDD(a_1, alpha=alpha)
            risk_2 = ConAbsDD(a_2, alpha=alpha)

        RC_i = (risk_1 - risk_2) / (2 * d_i) * w_[i, 0]
        RC.append(RC_i)

    RC = np.array(RC, ndmin=1)

    return RC
