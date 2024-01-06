""""""  #
"""
Copyright (c) 2020-2024, Dany Cajas
All rights reserved.
This work is licensed under BSD 3-Clause "New" or "Revised" License.
License available at https://github.com/dcajasn/Riskfolio-Lib/blob/master/LICENSE.txt
"""

import numpy as np
import pandas as pd
import cvxpy as cp
import scipy.stats as st
from scipy.linalg import sqrtm
import riskfolio.src.RiskFunctions as rk
import riskfolio.src.ParamsEstimation as pe
import riskfolio.src.AuxFunctions as af
import riskfolio.src.OwaWeights as owa
import riskfolio.external.cppfunctions as cf


__all__ = [
    "Portfolio",
]

rmeasures = [
    "MV",
    "KT",
    "MAD",
    "GMD",
    "MSV",
    "SKT",
    "FLPM",
    "SLPM",
    "CVaR",
    "TG",
    "EVaR",
    "RLVaR",
    "WR",
    "CVRG",
    "TGRG",
    "RG",
    "MDD",
    "ADD",
    "CDaR",
    "EDaR",
    "RLDaR",
    "UCI",
]


class Portfolio(object):
    r"""
    Class that creates a portfolio object with all properties needed to
    calculate optimal portfolios.

    Parameters
    ----------
    returns : DataFrame, optional
        A dataframe that containts the returns of the assets.
        The default is None.
    sht : bool, optional
        Indicate if the portfolio consider short positions (negative weights).
        The default is False.
    uppersht : float, optional
        Indicate the maximum value of the sum of absolute values of short
        positions (negative weights). The default is 0.2.
    upperlng : float, optional
        Indicate the maximum value of the sum of long positions (positive
        weights). When sht=True, the difference between upperlng and uppersht
        must be equal to the budget (upperlng - uppersht = budget)
        The default is 1.
    budget : float, optional
        Indicate the maximum value of the sum of long positions (positive
        weights) and short positions (negative weights). The default is 1.
    nea : int, optional
        Indicate the minimum number of effective assets (NEA) used in
        portfolio. This value is the inverse of Herfindahl-Hirschman index of
        portfolio's weights. The default is None.
    card : int, optional
        Indicate the maximum number of assets used in portfolio. It requires
        a solver that supports Mixed Integer Programs (MIP), see `Solvers <https://www.cvxpy.org/tutorial/advanced/index.html#solve-method-options>`_ for more details.
        This constraint is based on :cite:`a-YUE2014949`. The default is None.
    factors : DataFrame, optional
        A dataframe that containts the returns of the factors.
        The default is None.
    B : DataFrame, optional
        A dataframe that containts the loadings matrix.
        The default is None.
    alpha : float, optional
        Significance level of CVaR, EVaR, CDaR, EDaR and Tail Gini of losses. The default is 0.05.
    a_sim : float, optional
        Number of CVaRs used to approximate Tail Gini of losses. The default is 100.
    beta : float, optional
        Significance level of CVaR and Tail Gini of gains. If None it duplicates alpha value.
        The default is None.
    b_sim : float, optional
        Number of CVaRs used to approximate Tail Gini of gains. If None it duplicates a_sim value.
        The default is None.
    kappa : float, optional
        Deformation parameter of RLVaR and RLDaR, must be between 0 and 1. The default is 0.30.
    n_max_kurt : int, optional
        Maximum number of assets to use Kurtosis model based on semidefinte formulation. If number of
        assets is higher than n_max_kurt, it uses relaxed kurtosis model based on second order cone.
        The default is 50.
    kindbench : bool, optional
        True if the benchmark is a portfolio with detailed weights and False if
        the benchmark is an index. The default is True.
    allowTO : bool, optional
        Indicate if there is turnover constraints. The default is False.
    turnover : float, optional
        The maximum limit of turnover deviations. The default is 0.05.
    allowTE : bool, optional
        Indicate if there is tracking error constraints.. The default is False.
    TE : float, optional
        The maximum limit of tracking error deviations. The default is 0.05.
    benchindex : DataFrame, optional
        A dataframe that containts the returns of an index. If kindbench is
        False the tracking error constraints are calculated respect to this
        index. The default is None.
    benchweights : DataFrame, optional
        A dataframe that containts the weights of an index. The default is the
        equally weighted portfolio 1/N.
    ainequality : nd-array, optional
        The matrix :math:`A` of the linear constraint :math:`A \geq B`.
        The default is None.
    binequality : 1d-array, optional
        The matrix :math:`B` of the linear constraint :math:`A \geq B`.
        The default is None.
    b : 1d-array, optional
        The risk budgeting constraint vector. The default is None.
    network_sdp : nd-array, optional
        Connection matrix for semidefinite programming (SDP) network constraint.
        Users cannot use network_sdp and network_ip at the same time, when a value
        is assigned to network_sdp automatically network_ip becomes None.
        This constraint is based on :cite:`a-Cajas10` and :cite:`a-Cajas11`.
        The default is None.
    network_penalty : float, optional
        The weight of SDP network constraint when the risk measure is not 'MV'.
        This constraint is based on :cite:`a-Cajas10` and :cite:`a-Cajas11`.
        The default is 0.05.
    network_ip : nd-array, optional
        Connection matrix for Integer Programming (IP) network constraint. Users
        cannot use network_sdp and network_ip at the same time, when a value is
        assigned to network_ip automatically network_sdp becomes None.
        This constraint is based on :cite:`a-Cajas10` and :cite:`a-Cajas11`.
        The default is None.
    acentrality : nd-array, optional
        The matrix :math:`A_c` of the centrality constraint :math:`A_c = B_c`.
        This constraint is based on :cite:`a-Cajas10`.
        The default is None.
    bcentrality : nd-array, optional
        The matrix :math:`B_c` of the centrality constraint :math:`A_c = B_c`.
        This constraint is based on :cite:`a-Cajas10`.
        The default is None.
    lowerret : float, optional
        Constraint on min level of expected return. The default is None.
    upperdev : float, optional
        Constraint on max level of standard deviation. The default is None.
    upperkt : float, optional
        Constraint on max level of square root kurtosis. The default is None.
    uppermad : float, optional
        Constraint on max level of MAD. The default is None.
    uppergmd : float, optional
        Constraint on max level of GMD. The default is None.
    uppersdev : float, optional
        Constraint on max level of semi standard deviation. The default is None.
    upperskt : float, optional
        Constraint on max level of square root semi kurtosis. The default is None.
    upperflpm : float, optional
        Constraint on max level of first lower partial moment.
        The default is None.
    upperslpm : float, optional
        Constraint on max level of second lower partial moment.
        The default is None.
    upperCVaR : float, optional
        Constraint on max level of conditional value at risk (CVaR). The default is None.
    uppertg : float, optional
        Constraint on max level of Tail Gini. The default is None.
    upperEVaR : float, optional
        Constraint on max level of entropic value at risk (EVaR). The default is None.
    upperRLVaR : float, optional
        Constraint on max level of relativistic value at risk (RLVaR). The default is None.
    upperwr : float, optional
        Constraint on max level of worst realization. The default is None.
    upperrg : float, optional
        Constraint on max level of range. The default is None.
    uppercvrg : float, optional
        Constraint on max level of CVaR range. The default is None.
    uppertgrg : float, optional
        Constraint on max level of Tail Gini range. The default is None.
    uppermdd : float, optional
        Constraint on max level of maximum drawdown of uncompounded cumulative
        returns. The default is None.
    upperadd : float, optional
        Constraint on max level of average drawdown of uncompounded cumulative
        returns. The default is None.
    upperCDaR : float, optional
        Constraint on max level of conditional drawdown at risk (CDaR) of
        uncompounded cumulative returns. The default is None.
    upperEDaR : float, optional
        Constraint on max level of entropic drawdown at risk (EDaR) of
        uncompounded cumulative returns. The default is None.
    upperRLDaR : float, optional
        Constraint on max level of relativistic drawdown at risk (RLDaR) of
        uncompounded cumulative returns. The default is None.
    upperuci : float, optional
        Constraint on max level of ulcer index (UCI) of
        uncompounded cumulative returns. The default is None.
    """

    def __init__(
        self,
        returns=None,
        sht=False,
        uppersht=0.2,
        upperlng=1,
        budget=1,
        nea=None,
        card=None,
        factors=None,
        B=None,
        alpha=0.05,
        a_sim=100,
        beta=None,
        b_sim=None,
        kappa=0.30,
        n_max_kurt=50,
        kindbench=True,
        allowTO=False,
        turnover=0.05,
        allowTE=False,
        TE=0.05,
        benchindex=None,
        benchweights=None,
        ainequality=None,
        binequality=None,
        b=None,
        network_sdp=None,
        network_penalty=0.05,
        network_ip=None,
        acentrality=None,
        bcentrality=None,
        lowerret=None,
        upperdev=None,
        upperkt=None,
        uppermad=None,
        uppergmd=None,
        uppersdev=None,
        upperskt=None,
        upperflpm=None,
        upperslpm=None,
        upperCVaR=None,
        uppertg=None,
        upperEVaR=None,
        upperRLVaR=None,
        upperwr=None,
        uppercvrg=None,
        uppertgrg=None,
        upperrg=None,
        uppermdd=None,
        upperadd=None,
        upperCDaR=None,
        upperEDaR=None,
        upperRLDaR=None,
        upperuci=None,
    ):
        # Optimization Models Options

        self._returns = returns
        self.sht = sht
        self.uppersht = uppersht
        self.upperlng = upperlng
        self.budget = budget
        self.nea = nea
        self.card = card
        self._factors = factors
        self.alpha = alpha
        self.a_sim = a_sim
        self.beta = beta
        self.b_sim = b_sim
        self._kappa = kappa
        self.n_max_kurt = n_max_kurt
        self.kindbench = kindbench
        self.benchindex = benchindex
        self._benchweights = benchweights
        self._ainequality = ainequality
        self._binequality = binequality
        self._b = b
        self._network_sdp = network_sdp
        self.network_penalty = network_penalty
        self._network_ip = network_ip
        self._acentrality = acentrality
        self._bcentrality = bcentrality
        self.lowerret = lowerret
        self.upperdev = upperdev
        self.upperkt = upperkt
        self.uppermad = uppermad
        self.uppergmd = uppergmd
        self.uppersdev = uppersdev
        self.upperskt = upperskt
        self.upperflpm = upperflpm
        self.upperslpm = upperslpm
        self.upperCVaR = upperCVaR
        self.uppertg = uppertg
        self.upperEVaR = upperEVaR
        self.upperRLVaR = upperRLVaR
        self.upperwr = upperwr
        self.uppercvrg = uppercvrg
        self.uppertgrg = uppertgrg
        self.upperrg = upperrg
        self.uppermdd = uppermdd
        self.upperadd = upperadd
        self.upperCDaR = upperCDaR
        self.upperEDaR = upperEDaR
        self.upperRLDaR = upperRLDaR
        self.upperuci = upperuci

        self.allowTO = allowTO
        self.turnover = turnover
        self.allowTE = allowTE
        self.TE = TE

        # Inputs of Optimization Models

        self.mu = None
        self.cov = None
        self.kurt = None
        self.skurt = None
        self.L_2 = None
        self.S_2 = None
        self.mu_f = None
        self.cov_f = None
        self._B = None
        self.mu_fm = None
        self.cov_fm = None
        self.mu_bl = None
        self.cov_bl = None
        self.mu_bl_fm = None
        self.cov_bl_fm = None
        self.returns_fm = None
        self.z_EVaR = None
        self.z_EDaR = None
        self.z_RLVaR = None
        self.z_RLDaR = None

        # Inputs of Worst Case Optimization Models

        self.cov_l = None
        self.cov_u = None
        self.cov_mu = None
        self.cov_sigma = None
        self.d_mu = None
        self.k_mu = None
        self.k_sigma = None

        # Optimal portfolios

        self.optimal = None
        self.rp_optimal = None
        self.rrp_optimal = None
        self.wc_optimal = None
        self.limits = None
        self.frontier = None

        # Solver params

        self.solvers = ["CLARABEL", "ECOS", "SCS", "OSQP", "CVXOPT"]
        self.sol_params = {
            # 'ECOS': {"max_iters": 500, "abstol": 1e-8},
            # 'SCS': {"max_iters": 2500, "eps": 1e-5},
            # 'OSQP': {"max_iter": 10000, "eps_abs": 1e-8},
            # 'CVXOPT': {"max_iters": 500, "abstol": 1e-8},
        }

    @property
    def returns(self):
        if self._returns is not None and isinstance(self._returns, pd.DataFrame):
            return self._returns
        else:
            raise NameError("returns must be a DataFrame")

    @returns.setter
    def returns(self, value):
        if value is not None and isinstance(value, pd.DataFrame):
            self._returns = value
        else:
            raise NameError("returns must be a DataFrame")

    @property
    def assetslist(self):
        if self._returns is not None and isinstance(self._returns, pd.DataFrame):
            return self._returns.columns.tolist()
        elif self._returns is None:
            return None

    @property
    def numassets(self):
        if self._returns is not None and isinstance(self._returns, pd.DataFrame):
            return self._returns.shape[1]

    @property
    def factors(self):
        return self._factors

    @factors.setter
    def factors(self, value):
        a = value
        if a is not None and isinstance(a, pd.DataFrame):
            if self.returns.index.equals(a.index):
                self._factors = a
        else:
            raise NameError("factors must be a DataFrame")

    @property
    def factorslist(self):
        if self._factors is not None and isinstance(self._factors, pd.DataFrame):
            return self._factors.columns.tolist()
        elif self._factors is None:
            return None

    @property
    def B(self):
        return self._B

    @B.setter
    def B(self, value):
        a = value
        if a is not None and isinstance(a, pd.DataFrame):
            self._B = a
        else:
            raise NameError("loadings matrix must be a DataFrame")

    @property
    def benchweights(self):
        n = self.numassets
        if self._benchweights is not None:
            if self._benchweights.shape[0] == n and self._benchweights.shape[1] == 1:
                a = self._benchweights
            else:
                raise NameError("Weights must have a size of shape (n_assets,1)")
        else:
            a = np.array(np.ones((n, 1)) / n)
        return a

    @benchweights.setter
    def benchweights(self, value):
        a = value
        n = self.numassets
        if a is not None:
            if a.shape[0] == n and a.shape[1] == 1:
                a = a
            else:
                raise NameError("Weights must have a size of shape (n_assets,1)")
        else:
            a = np.array(np.ones((n, 1)) / n)
        self._benchweights = a

    @property
    def ainequality(self):
        a = self._ainequality
        if a is not None:
            if a.shape[1] == self.numassets:
                a = a
            else:
                raise NameError(
                    "The array ainequality must have the same number of columns than assets' number"
                )
        return a

    @ainequality.setter
    def ainequality(self, value):
        a = value
        if a is not None:
            if a.shape[1] == self.numassets:
                a = a
            else:
                raise NameError(
                    "The matrix ainequality must have the same number of columns than assets' number"
                )
        self._ainequality = a

    @property
    def binequality(self):
        a = self._binequality
        if a is not None:
            if a.shape[1] == 1:
                a = a
            else:
                raise NameError("The matrix binequality must have one column")
        return a

    @binequality.setter
    def binequality(self, value):
        a = value
        if a is not None:
            if a.shape[1] == 1:
                a = a
            else:
                raise NameError("The matrix binequality must have one column")
        self._binequality = a

    @property
    def b(self):
        a = self._b
        if a is not None:
            if a.shape[0] == self.numassets and a.shape[1] == 1:
                pass
            elif a.shape[0] == 1 and a.shape[1] == self.numassets:
                a = a.T
            else:
                raise NameError(
                    "The vector of risk contribution constraints must have a size equal than the assets' number"
                )
        return a

    @b.setter
    def b(self, value):
        a = value
        if a is not None:
            if a.shape[0] == self.numassets and a.shape[1] == 1:
                pass
            elif a.shape[0] == 1 and a.shape[1] == self.numassets:
                a = a.T
            else:
                raise NameError(
                    "The vector of risk contribution constraints must have a size equal than the assets' number"
                )
        self._b = a

    @property
    def network_sdp(self):
        a = self._network_sdp
        n = self.numassets
        if self._network_sdp is not None:
            if self._network_sdp.shape[0] == n and self._network_sdp.shape[1] == n:
                a = self._network_sdp
            else:
                raise NameError(
                    "Connection matrix network_sdp must have a size of shape (n_assets,n_assets)"
                )
        return a

    @network_sdp.setter
    def network_sdp(self, value):
        a = value
        n = self.numassets
        if a is not None:
            if a.shape[0] == n and a.shape[1] == n:
                a = a
            else:
                raise NameError(
                    "Connection matrix network_sdp must have a size of shape (n_assets,n_assets)"
                )
        self._network_sdp = a
        self._network_ip = None

    @property
    def network_ip(self):
        a = self._network_ip
        n = self.numassets
        if a is not None:
            if a.shape[1] == n:
                a = a
            else:
                raise NameError(
                    "Connection matrix network_ip must have a n_assets columns"
                )
        return a

    @network_ip.setter
    def network_ip(self, value):
        a = value
        n = self.numassets
        if a is not None:
            if a.shape[1] == n:
                a = a
            else:
                raise NameError(
                    "Connection matrix network_ip must have a n_assets columns"
                )
        self._network_ip = a
        self._network_sdp = None

    @property
    def acentrality(self):
        a = self._acentrality
        if a is not None:
            if a.shape[1] == self.numassets:
                a = a
            else:
                raise NameError(
                    "The array ainequality must have the same number of columns than assets' number"
                )
        return a

    @acentrality.setter
    def acentrality(self, value):
        a = value
        if a is not None:
            if a.shape[1] == self.numassets:
                a = a
            else:
                raise NameError(
                    "The matrix ainequality must have the same number of columns than assets' number"
                )
        self._acentrality = a

    @property
    def bcentrality(self):
        a = self._bcentrality
        if a is not None:
            if a.shape[1] == 1:
                a = a
            else:
                raise NameError("The matrix binequality must have one column")
        return a

    @bcentrality.setter
    def bcentrality(self, value):
        a = value
        if a is not None:
            if a.shape[1] == 1:
                a = a
            else:
                raise NameError("The matrix binequality must have one column")
        self._bcentrality = a

    @property
    def kappa(self):
        return self._kappa

    @kappa.setter
    def kappa(self, value):
        a = value
        if a >= 1:
            print(
                "kappa must be between 0 and 1, values higher or equal to 1 are setting to 0.99"
            )
            self._kappa = 0.99
        elif a <= 0:
            print(
                "kappa must be between 0 and 1, values lower or equal to 0 are setting to 0.01"
            )
            self._kappa = 0.01
        else:
            self._kappa = a

    def assets_stats(
        self, method_mu="hist", method_cov="hist", method_kurt=None, d=0.94, **kwargs
    ):
        r"""
        Calculate the inputs that will be used by the optimization method when
        we select the input model='Classic'.

        Parameters
        ----------
        method_mu : str, optional
            The method used to estimate the expected returns.
            The default value is 'hist'. Possible values are:

            - 'hist': use historical estimates.
            - 'ewma1'': use ewma with adjust=True, see `EWM <https://pandas.pydata.org/pandas-docs/stable/user_guide/computation.html#exponentially-weighted-windows>`_ for more details.
            - 'ewma2': use ewma with adjust=False, see `EWM <https://pandas.pydata.org/pandas-docs/stable/user_guide/computation.html#exponentially-weighted-windows>`_ for more details.
            - 'JS': James-Stein estimator. For more information see :cite:`a-Meucci2005` and :cite:`a-Feng2016`.
            - 'BS': Bayes-Stein estimator. For more information see :cite:`a-Jorion1986`.
            - 'BOP': BOP estimator. For more information see :cite:`a-Bodnar2019`.

        method_cov : str, optional
            The method used to estimate the covariance matrix.
            The default is 'hist'. Possible values are:

            - 'hist': use historical estimates.
            - 'ewma1'': use ewma with adjust=True, see `EWM <https://pandas.pydata.org/pandas-docs/stable/user_guide/computation.html#exponentially-weighted-windows>`_ for more details.
            - 'ewma2': use ewma with adjust=False, see `EWM <https://pandas.pydata.org/pandas-docs/stable/user_guide/computation.html#exponentially-weighted-windows>`_ for more details.
            - 'ledoit': use the Ledoit and Wolf Shrinkage method.
            - 'oas': use the Oracle Approximation Shrinkage method.
            - 'shrunk': use the basic Shrunk Covariance method.
            - 'gl': use the basic Graphical Lasso Covariance method.
            - 'jlogo': use the j-LoGo Covariance method. For more information see: :cite:`a-jLogo`.
            - 'fixed': denoise using fixed method. For more information see chapter 2 of :cite:`a-MLforAM`.
            - 'spectral': denoise using spectral method. For more information see chapter 2 of :cite:`a-MLforAM`.
            - 'shrink': denoise using shrink method. For more information see chapter 2 of :cite:`a-MLforAM`.
            - 'gerber1': use the Gerber statistic 1. For more information see: :cite:`a-Gerber2021`.
            - 'gerber2': use the Gerber statistic 2. For more information see: :cite:`a-Gerber2021`.

        method_kurt : str, optional
            The method used to estimate the kurtosis square matrix:
            The default is None. Possible values are:

            - None: do not calculate kurtosis square matrix.
            - 'hist': use historical estimates. For more information see :cite:`a-Cajas4`.
            - 'semi': use semi cokurtosis square matrix. For more information see :cite:`a-Cajas4`.
            - 'fixed': denoise using fixed method. For more information see chapter 2 of :cite:`a-MLforAM`.
            - 'spectral': denoise using spectral method. For more information see chapter 2 of :cite:`a-MLforAM`.
            - 'shrink': denoise using shrink method. For more information see chapter 2 of :cite:`a-MLforAM`.

        **kwargs : dict
            All additional parameters of mean_vector and covar_matrix functions.

        See Also
        --------
        riskfolio.src.ParamsEstimation.mean_vector
        riskfolio.src.ParamsEstimation.covar_matrix
        riskfolio.src.ParamsEstimation.cokurt_matrix

        """

        self.mu = pe.mean_vector(self.returns, method=method_mu, d=d)
        self.cov = pe.covar_matrix(self.returns, method=method_cov, **kwargs)
        value = af.is_pos_def(self.cov, threshold=1e-8)
        for i in range(5):
            if value == False:
                try:
                    self.cov = af.cov_fix(self.cov, method="clipped", threshold=1e-5)
                    value = af.is_pos_def(self.cov, threshold=1e-8)
                except:
                    break
            else:
                break

        if value == False:
            print("You must convert self.cov to a positive definite matrix")

        if method_kurt not in [None, "semi"]:
            T, N = self.returns.shape
            self.L_2 = cf.duplication_elimination_matrix(N)
            self.S_2 = cf.duplication_summation_matrix(N)
            self.kurt = pe.cokurt_matrix(self.returns, method=method_kurt, **kwargs)
            value = af.is_pos_def(self.kurt, threshold=1e-8)
            for i in range(5):
                if value == False:
                    try:
                        self.kurt = af.cov_fix(
                            self.kurt, method="clipped", threshold=1e-5
                        )
                        value = af.is_pos_def(self.kurt, threshold=1e-8)
                    except:
                        break
                else:
                    break

            if value == False:
                print("You must convert self.kurt to a positive definite matrix")

            self.skurt = pe.cokurt_matrix(self.returns, method="semi")
            value = af.is_pos_def(self.skurt, threshold=1e-8)
            for i in range(5):
                if value == False:
                    try:
                        self.skurt = af.cov_fix(
                            self.skurt, method="clipped", threshold=1e-5
                        )
                        value = af.is_pos_def(self.skurt, threshold=1e-8)
                    except:
                        break
                else:
                    break

            if value == False:
                print("You must convert self.skurt to a positive definite matrix")

        else:
            self.kurt = None
            self.skurt = None
            self.L_2 = None
            self.S_2 = None

    def blacklitterman_stats(
        self,
        P,
        Q,
        rf=0,
        w=None,
        delta=None,
        eq=True,
        method_mu="hist",
        method_cov="hist",
        **kwargs,
    ):
        r"""
        Calculate the inputs that will be used by the optimization method when
        we select the input model='BL'.

        Parameters
        ----------
        P : DataFrame of shape (n_views, n_assets)
            Analyst's views matrix, can be relative or absolute.
        Q: DataFrame of shape (n_views, 1)
            Expected returns of analyst's views.
        delta: float
            Risk aversion factor. The default value is 1.
        rf: scalar, optional
            Risk free rate. The default is 0.
        w : DataFrame of shape (n_assets, 1)
            Weights matrix, where n_assets is the number of assets.
            The default is None.
        eq: bool, optional
            Indicates if use equilibrium or historical excess returns.
            The default is True.
        method_mu : str, optional
            The method used to estimate the expected returns.
            The default value is 'hist'. Possible values are:

            - 'hist': use historical estimates.
            - 'ewma1'': use ewma with adjust=True, see `EWM <https://pandas.pydata.org/pandas-docs/stable/user_guide/computation.html#exponentially-weighted-windows>`_ for more details.
            - 'ewma2': use ewma with adjust=False, see `EWM <https://pandas.pydata.org/pandas-docs/stable/user_guide/computation.html#exponentially-weighted-windows>`_ for more details.
            - 'JS': James-Stein estimator. For more information see :cite:`a-Meucci2005` and :cite:`a-Feng2016`.
            - 'BS': Bayes-Stein estimator. For more information see :cite:`a-Jorion1986`.
            - 'BOP': BOP estimator. For more information see :cite:`a-Bodnar2019`.
        method_cov : str, optional
            The method used to estimate the covariance matrix.
            The default is 'hist'. Possible values are:

            - 'hist': use historical estimates.
            - 'ewma1'': use ewma with adjust=True, see `EWM <https://pandas.pydata.org/pandas-docs/stable/user_guide/computation.html#exponentially-weighted-windows>`_ for more details.
            - 'ewma2': use ewma with adjust=False, see `EWM <https://pandas.pydata.org/pandas-docs/stable/user_guide/computation.html#exponentially-weighted-windows>`_ for more details.
            - 'ledoit': use the Ledoit and Wolf Shrinkage method.
            - 'oas': use the Oracle Approximation Shrinkage method.
            - 'shrunk': use the basic Shrunk Covariance method.
            - 'gl': use the basic Graphical Lasso Covariance method.
            - 'jlogo': use the j-LoGo Covariance method. For more information see: :cite:`a-jLogo`.
            - 'fixed': denoise using fixed method. For more information see chapter 2 of :cite:`a-MLforAM`.
            - 'spectral': denoise using spectral method. For more information see chapter 2 of :cite:`a-MLforAM`.
            - 'shrink': denoise using shrink method. For more information see chapter 2 of :cite:`a-MLforAM`.
            - 'gerber1': use the Gerber statistic 1. For more information see: :cite:`a-Gerber2021`.
            - 'gerber2': use the Gerber statistic 2. For more information see: :cite:`a-Gerber2021`.
        **kwargs : dict
            Other variables related to the covariance estimation.

        See Also
        --------
        riskfolio.src.ParamsEstimation.black_litterman

        """
        X = self.returns
        if w is None:
            w = np.array(self.benchweights, ndmin=2)

        if delta is None:
            a = np.array(self.mu, ndmin=2) @ np.array(w, ndmin=2)
            delta = (a - rf) / (
                np.array(w, ndmin=2).T
                @ np.array(self.cov, ndmin=2)
                @ np.array(w, ndmin=2)
            )
            delta = delta.item()

        mu, cov, w = pe.black_litterman(
            X=X,
            w=w,
            P=P,
            Q=Q,
            delta=delta,
            rf=rf,
            eq=eq,
            method_mu=method_mu,
            method_cov=method_cov,
            **kwargs,
        )
        self.mu_bl = mu
        self.cov_bl = cov

        value = af.is_pos_def(self.cov_bl, threshold=1e-8)
        for i in range(5):
            if value == False:
                try:
                    self.cov_bl = af.cov_fix(
                        self.cov_bl, method="clipped", threshold=1e-5
                    )
                    value = af.is_pos_def(self.cov_bl, threshold=1e-8)
                except:
                    break
            else:
                break

        if value == False:
            print("You must convert self.cov_bl to a positive definite matrix")

    def factors_stats(
        self,
        method_mu="hist",
        method_cov="hist",
        d=0.94,
        B=None,
        dict_cov={},
        dict_risk={},
    ):
        r"""
        Calculate the inputs that will be used by the optimization method when
        we select the input model='FM'.

        Parameters
        ----------
        method_mu : str, optional
            The method used to estimate the expected returns.
            The default value is 'hist'. Possible values are:

            - 'hist': use historical estimates.
            - 'ewma1'': use ewma with adjust=True, see `EWM <https://pandas.pydata.org/pandas-docs/stable/user_guide/computation.html#exponentially-weighted-windows>`_ for more details.
            - 'ewma2': use ewma with adjust=False, see `EWM <https://pandas.pydata.org/pandas-docs/stable/user_guide/computation.html#exponentially-weighted-windows>`_ for more details.
            - 'JS': James-Stein estimator. For more information see :cite:`a-Meucci2005` and :cite:`a-Feng2016`.
            - 'BS': Bayes-Stein estimator. For more information see :cite:`a-Jorion1986`.
            - 'BOP': BOP estimator. For more information see :cite:`a-Bodnar2019`.
        method_cov : str, optional
            The method used to estimate the covariance matrix.
            The default is 'hist'. Possible values are:

            - 'hist': use historical estimates.
            - 'ewma1'': use ewma with adjust=True, see `EWM <https://pandas.pydata.org/pandas-docs/stable/user_guide/computation.html#exponentially-weighted-windows>`_ for more details.
            - 'ewma2': use ewma with adjust=False, see `EWM <https://pandas.pydata.org/pandas-docs/stable/user_guide/computation.html#exponentially-weighted-windows>`_ for more details.
            - 'ledoit': use the Ledoit and Wolf Shrinkage method.
            - 'oas': use the Oracle Approximation Shrinkage method.
            - 'shrunk': use the basic Shrunk Covariance method.
            - 'gl': use the basic Graphical Lasso Covariance method.
            - 'jlogo': use the j-LoGo Covariance method. For more information see: :cite:`a-jLogo`.
            - 'fixed': denoise using fixed method. For more information see chapter 2 of :cite:`a-MLforAM`.
            - 'spectral': denoise using spectral method. For more information see chapter 2 of :cite:`a-MLforAM`.
            - 'shrink': denoise using shrink method. For more information see chapter 2 of :cite:`a-MLforAM`.
            - 'gerber1': use the Gerber statistic 1. For more information see: :cite:`a-Gerber2021`.
            - 'gerber2': use the Gerber statistic 2. For more information see: :cite:`a-Gerber2021`.
        dict_cov : dict
            Other variables related to the covariance estimation.
        dict_risk : dict
            Other variables related of risk_factors function.

        See Also
        --------
        riskfolio.src.ParamsEstimation.forward_regression
        riskfolio.src.ParamsEstimation.backward_regression
        riskfolio.src.ParamsEstimation.loadings_matrix
        riskfolio.src.ParamsEstimation.risk_factors

        """
        X = self.factors
        Y = self.returns

        mu_f = pe.mean_vector(self.returns, method=method_mu, d=d)
        cov_f = pe.covar_matrix(self.returns, method=method_cov, d=d, **dict_cov)

        self.mu_f = mu_f
        self.cov_f = cov_f

        value = af.is_pos_def(self.cov_f, threshold=1e-8)
        if value == False:
            try:
                self.cov = af.cov_fix(self.cov, method="clipped", threshold=1e-5)
                value = af.is_pos_def(self.cov, threshold=1e-8)
                if value == False:
                    print("You must convert self.cov to a positive definite matrix")
            except:
                print("You must convert self.cov to a positive definite matrix")

        mu, cov, returns = pe.risk_factors(
            X, Y, B=B, method_mu=method_mu, method_cov=method_cov, **dict_risk
        )

        self.mu_fm = mu
        self.cov_fm = cov
        self.returns_fm = returns

        value = af.is_pos_def(self.cov_fm, threshold=1e-8)
        for i in range(5):
            if value == False:
                try:
                    self.cov_fm = af.cov_fix(
                        self.cov_fm, method="clipped", threshold=1e-5
                    )
                    value = af.is_pos_def(self.cov_fm, threshold=1e-8)
                except:
                    break
            else:
                break

        if value == False:
            print("You must convert self.cov_fm to a positive definite matrix")

    def blfactors_stats(
        self,
        flavor="BLB",
        B=None,
        P=None,
        Q=None,
        P_f=None,
        Q_f=None,
        rf=0,
        w=None,
        delta=None,
        eq=True,
        const=False,
        diag=False,
        method_mu="hist",
        method_cov="hist",
        kwargs_1=None,
        kwargs_2=None,
    ):
        r"""
        Calculate the inputs that will be used by the optimization method when
        we select the input model='BL'.

        Parameters
        ----------
        flavor : str
            Model used, can be 'BLB' for Black Litterman Bayesian or 'ABL' for
            Augmented Black Litterman. The default value is 'BLB'.
        B : DataFrame of shape (n_assets, n_features)
            Loadings matrix. The default value is None.
        P : DataFrame of shape (n_views, n_assets)
            Analyst's views matrix, can be relative or absolute.
        Q: DataFrame of shape (n_views, 1)
            Expected returns of analyst's views.
        P_f : DataFrame of shape (n_views, n_assets)
            Analyst's factors views matrix, can be relative or absolute.
        Q_f: DataFrame of shape (n_views, 1)
            Expected returns of analyst's factors views.
        delta: float
            Risk aversion factor. The default value is 1.
        rf: scalar, optional
            Risk free rate. The default is 0.
        w : DataFrame of shape (n_assets, 1)
            Weights matrix, where n_assets is the number of assets.
            The default is None.
        eq: bool, optional
            Indicates if use equilibrium or historical excess returns.
            The default is True.
        const : bool, optional
            Indicate if the loadings matrix has a constant.
            The default is False.
        diag : bool, optional
            Indicate if we use the diagonal matrix to calculate covariance matrix
            of factor model, only useful when we work with a factor model based on
            a regresion model (only equity portfolio).
            The default is False.
        method_mu : str, optional
            The method used to estimate the expected returns.
            The default value is 'hist'. Possible values are:

            - 'hist': use historical estimates.
            - 'ewma1'': use ewma with adjust=True, see `EWM <https://pandas.pydata.org/pandas-docs/stable/user_guide/computation.html#exponentially-weighted-windows>`_ for more details.
            - 'ewma2': use ewma with adjust=False, see `EWM <https://pandas.pydata.org/pandas-docs/stable/user_guide/computation.html#exponentially-weighted-windows>`_ for more details.
            - 'JS': James-Stein estimator. For more information see :cite:`a-Meucci2005` and :cite:`a-Feng2016`.
            - 'BS': Bayes-Stein estimator. For more information see :cite:`a-Jorion1986`.
            - 'BOP': BOP estimator. For more information see :cite:`a-Bodnar2019`.
        method_cov : str, optional
            The method used to estimate the covariance matrix:
            The default is 'hist'. Possible values are:

            - 'hist': use historical estimates.
            - 'ewma1'': use ewma with adjust=True, see `EWM <https://pandas.pydata.org/pandas-docs/stable/user_guide/computation.html#exponentially-weighted-windows>`_ for more details.
            - 'ewma2': use ewma with adjust=False, see `EWM <https://pandas.pydata.org/pandas-docs/stable/user_guide/computation.html#exponentially-weighted-windows>`_ for more details.
            - 'ledoit': use the Ledoit and Wolf Shrinkage method.
            - 'oas': use the Oracle Approximation Shrinkage method.
            - 'shrunk': use the basic Shrunk Covariance method.
            - 'gl': use the basic Graphical Lasso Covariance method.
            - 'jlogo': use the j-LoGo Covariance method. For more information see: :cite:`a-jLogo`.
            - 'fixed': denoise using fixed method. For more information see chapter 2 of :cite:`a-MLforAM`.
            - 'spectral': denoise using spectral method. For more information see chapter 2 of :cite:`a-MLforAM`.
            - 'shrink': denoise using shrink method. For more information see chapter 2 of :cite:`a-MLforAM`.
            - 'gerber1': use the Gerber statistic 1. For more information see: :cite:`a-Gerber2021`.
            - 'gerber2': use the Gerber statistic 2. For more information see: :cite:`a-Gerber2021`.
        kwargs_1 : dict
            Other variables related to the loadings matrix estimation.
        kwargs_2 : dict
            Other variables related to the factors Black Litterman model selected.

        See Also
        --------
        riskfolio.src.ParamsEstimation.augmented_black_litterman
        riskfolio.src.ParamsEstimation.black_litterman_bayesian

        """
        X = self.returns
        F = self.factors

        if w is None:
            w = np.array(self.benchweights, ndmin=2)

        if delta is None:
            a = np.array(self.mu, ndmin=2) @ np.array(w, ndmin=2)
            delta = (a - rf) / (
                np.array(w, ndmin=2).T
                @ np.array(self.cov, ndmin=2)
                @ np.array(w, ndmin=2)
            )
            delta = delta.item()

        if B is None:
            if self.B is None:
                self.B = pe.loadings_matrix(X=F, Y=X, **kwargs_1)
                const = True
            elif self.B is not None:
                pass
        elif B is not None:
            self.B = B

        if flavor == "BLB":
            if isinstance(kwargs_1, dict):
                mu, cov, w = pe.black_litterman_bayesian(
                    X=X,
                    F=F,
                    B=self.B,
                    P_f=P_f,
                    Q_f=Q_f,
                    delta=delta,
                    rf=rf,
                    eq=eq,
                    const=const,
                    diag=diag,
                    method_mu=method_mu,
                    method_cov=method_cov,
                    **kwargs_2,
                )
            else:
                mu, cov, w = pe.black_litterman_bayesian(
                    X=X,
                    F=F,
                    B=self.B,
                    P_f=P_f,
                    Q_f=Q_f,
                    delta=delta,
                    rf=rf,
                    eq=eq,
                    const=const,
                    diag=diag,
                    method_mu=method_mu,
                    method_cov=method_cov,
                )

        elif flavor == "ABL":
            if isinstance(kwargs_1, dict):
                mu, cov, w = pe.augmented_black_litterman(
                    X=X,
                    w=w,
                    F=F,
                    B=self.B,
                    P=P,
                    Q=Q,
                    P_f=P_f,
                    Q_f=Q_f,
                    delta=delta,
                    rf=rf,
                    eq=eq,
                    const=const,
                    method_mu=method_mu,
                    method_cov=method_cov,
                    **kwargs_2,
                )
            else:
                mu, cov, w = pe.augmented_black_litterman(
                    X=X,
                    w=w,
                    F=F,
                    B=self.B,
                    P=P,
                    Q=Q,
                    P_f=P_f,
                    Q_f=Q_f,
                    delta=delta,
                    rf=rf,
                    eq=eq,
                    const=const,
                    method_mu=method_mu,
                    method_cov=method_cov,
                )

        self.mu_bl_fm = mu
        self.cov_bl_fm = cov

        value = af.is_pos_def(self.cov_bl_fm, threshold=1e-8)
        for i in range(5):
            if value == False:
                try:
                    self.cov_bl_fm = af.cov_fix(
                        self.cov_bl_fm, method="clipped", threshold=1e-5
                    )
                    value = af.is_pos_def(self.cov_bl_fm, threshold=1e-8)
                except:
                    break
            else:
                break

        if value == False:
            print("You must convert self.cov_bl_fm to a positive definite matrix")

    def wc_stats(
        self,
        box="s",
        ellip="s",
        q=0.05,
        n_sim=3000,
        window=3,
        dmu=0.1,
        dcov=0.1,
        method_k_mu="normal",
        method_k_sigma="normal",
        seed=0,
    ):
        r"""
        Calculate the inputs that will be used by the wc_optimization method.

        Parameters
        ----------
        box : string
            The method used to estimate the box uncertainty sets. The default is 's'. Possible values are:

            - 's': stationary bootstrapping method.
            - 'c': circular bootstrapping method.
            - 'm': moving bootstrapping method.
            - 'n': assuming normal returns to calculate confidence levels.
            - 'd': delta method, this method increase and decrease by a percentage.

        ellip : string
            The method used to estimate the elliptical uncertainty sets. The default is 's'. Possible values are:

            - 's': stationary bootstrapping method.
            - 'c': circular bootstrapping method.
            - 'm': moving bootstrapping method.
            - 'n': assuming normal returns to calculate confidence levels.

        q : scalar
            Significance level of the selected bootstrapping method.
            The default is 0.05.
        n_sim : scalar
            Number of simulations of the bootstrapping method.
            The default is 3000.
        window:
            Block size of the bootstrapping method. Must be greather than 1
            and lower than the n_samples - n_features + 1
            The default is 3.
        dmu : scalar
            Percentage used by delta method to increase and decrease the mean vector in box constraints.
            The default is 0.1.
        dcov : scalar
            Percentage used by delta method to increase and decrease the covariance matrix in box constraints.
            The default is 0.1.
        method_k_mu : string
            Method used to calculate the elliptical uncertainty set of the mean.
            The default is 'normal'. Possible values are:

                - 'normal': assumes normal distribution of returns. Uses a Chi2 distribution. :cite:`a-Fabozzi`.
                - 'general': for any possible distribution of returns. Uses the ratio √((1-q)/q).  :cite:`a-Fabozzi`.

        method_k_sigma : string
            Method used to calculate the elliptical uncertainty set of the covariance matrix.
            The default is 'normal'. Possible values are:

                - 'normal': assumes normal distribution of returns. Uses a Chi2 distribution. :cite:`a-Fabozzi`.
                - 'general': for any possible distribution of returns. Uses the ratio √((1-q)/q).  :cite:`a-Fabozzi`,

        seed: int
            Seed used to generate the boostrapping sample. The defailt is 0.

        See Also
        --------
        riskfolio.src.ParamsEstimation.bootstrapping

        """

        if box not in list("scmdn"):
            raise ValueError("box only can be 's', 'c', 'm', 'd' or 'n'")
        if ellip not in list("scmn"):
            raise ValueError("box only can be 's', 'c', 'm' or 'n'")

        X = self.returns
        cols = X.columns.tolist()
        cols_2 = [i + "-" + j for i in cols for j in cols]
        (T, N) = X.shape
        mu = X.mean().to_frame().T
        cov = X.cov()

        if box == "s":
            mu_l, mu_u, cov_l, cov_u, _, _ = pe.bootstrapping(
                X, kind="stationary", q=q, n_sim=n_sim, window=window, seed=seed
            )
            d_mu = (mu_u - mu_l) / 2
        elif box == "c":
            mu_l, mu_u, cov_l, cov_u, _, _ = pe.bootstrapping(
                X, kind="circular", q=q, n_sim=n_sim, window=window, seed=seed
            )
            d_mu = (mu_u - mu_l) / 2
        elif box == "m":
            mu_l, mu_u, cov_l, cov_u, _, _ = pe.bootstrapping(
                X, kind="moving", q=q, n_sim=n_sim, window=window, seed=seed
            )
            d_mu = (mu_u - mu_l) / 2
        elif box == "n":
            # Defining confidence level of mean vector assuming normal returns
            d_mu = st.norm.ppf(1 - q / 2) * np.sqrt(np.diag(cov) / T).reshape(1, -1)
            d_mu = pd.DataFrame(d_mu, index=[0], columns=cols)

            # Defining confidence level of covariance matrix assuming normal returns
            rs = np.random.RandomState(seed=seed)
            A = st.wishart.rvs(T, cov / T, size=10000, random_state=rs)
            cov_l = np.percentile(A, q=q / 2, axis=0)
            cov_u = np.percentile(A, q=1 - q / 2, axis=0)

            cov_l = pd.DataFrame(cov_l, index=cols, columns=cols)
            cov_u = pd.DataFrame(cov_u, index=cols, columns=cols)

            if af.is_pos_def(cov_l) == False:
                cov_l = af.cov_fix(cov_l, method="clipped", threshold=1e-3)

            if af.is_pos_def(cov_u) == False:
                cov_u = af.cov_fix(cov_u, method="clipped", threshold=1e-3)

        elif box == "d":
            d_mu = dmu * np.abs(mu)
            cov_l = cov - dcov * np.abs(cov)
            cov_u = cov + dcov * np.abs(cov)
            d_mu = pd.DataFrame(d_mu, index=[0], columns=cols)
            cov_l = pd.DataFrame(cov_l, index=cols, columns=cols)
            cov_u = pd.DataFrame(cov_u, index=cols, columns=cols)

        if ellip == "s":
            _, _, _, _, cov_mu, cov_sigma = pe.bootstrapping(
                X, kind="stationary", q=q, n_sim=n_sim, window=window, seed=seed
            )
        elif ellip == "c":
            _, _, _, _, cov_mu, cov_sigma = pe.bootstrapping(
                X, kind="circular", q=q, n_sim=n_sim, window=window, seed=seed
            )
        elif ellip == "m":
            _, _, _, _, cov_mu, cov_sigma = pe.bootstrapping(
                X, kind="moving", q=q, n_sim=n_sim, window=window, seed=seed
            )
        elif ellip == "n":
            # Covariance of mean returns
            cov_mu = cov / T
            cov_mu = np.diag(np.diag(cov_mu))
            cov_mu = pd.DataFrame(cov_mu, index=cols, columns=cols)
            # Covariance of covariance matrix
            K = cf.commutation_matrix(N, N)
            I = np.identity(N * N)
            cov_sigma = T * (I + K) @ np.kron(cov_mu, cov_mu)
            cov_sigma = np.diag(np.diag(cov_sigma))
            cov_sigma = pd.DataFrame(cov_sigma, index=cols_2, columns=cols_2)

        self.cov_l = cov_l
        self.cov_u = cov_u
        self.cov_mu = cov_mu
        self.cov_sigma = cov_sigma
        self.d_mu = d_mu

        if method_k_mu == "normal":
            self.k_mu = st.chi2.ppf(1 - q, df=N) ** 0.5
        elif method_k_mu == "general":
            self.k_mu = np.sqrt((1 - q) / q)
        else:
            raise ValueError(
                "The only available values of parameter method_k_mu are 'normal' and 'general'"
            )

        if method_k_sigma == "normal":
            self.k_sigma = st.chi2.ppf(1 - q, df=N * N) ** 0.5
        elif method_k_sigma == "general":
            self.k_sigma = np.sqrt((1 - q) / q)
        else:
            raise ValueError(
                "The only available values of parameter method_k_mu are 'normal' and 'general'"
            )

    def optimization(
        self, model="Classic", rm="MV", obj="Sharpe", kelly=False, rf=0, l=2, hist=True
    ):
        r"""
        This method that calculates the optimal portfolio according to the
        optimization model selected by the user. The general problem that
        solves is:

        .. math::
            \begin{align}
            &\underset{w}{\text{optimize}} & & F(w)\\
            &\text{s. t.} & & Aw \geq B\\
            & & & \phi_{i}(w) \leq c_{i}\\
            \end{align}

        Where:

        :math:`F(w)` is the objective function.

        :math:`Aw \geq B` is a set of linear constraints.

        :math:`\phi_{i}(w) \leq c_{i}` are constraints on maximum values of
        several risk measures.

        Parameters
        ----------
        model : str can be {'Classic', 'BL', 'FM' or 'BLFM'}
            The model used for optimize the portfolio.
            The default is 'Classic'. Possible values are:

            - 'Classic': use estimates of expected return vector and covariance matrix that depends on historical data.
            - 'BL': use estimates of expected return vector and covariance matrix based on the Black Litterman model.
            - 'FM': use estimates of expected return vector and covariance matrix based on a Risk Factor model specified by the user.
            - 'BLFM': use estimates of expected return vector and covariance matrix based on Black Litterman applied to a Risk Factor model specified by the user.

        rm : str, optional
            The risk measure used to optimize the portfolio.
            The default is 'MV'. Possible values are:

            - 'MV': Standard Deviation.
            - 'KT': Square Root of Kurtosis.
            - 'MAD': Mean Absolute Deviation.
            - 'GMD': Gini Mean Difference.
            - 'MSV': Semi Standard Deviation.
            - 'SKT': Square Root of Semi Kurtosis.
            - 'FLPM': First Lower Partial Moment (Omega Ratio).
            - 'SLPM': Second Lower Partial Moment (Sortino Ratio).
            - 'CVaR': Conditional Value at Risk.
            - 'TG': Tail Gini.
            - 'EVaR': Entropic Value at Risk.
            - 'RLVaR': Relativistic Value at Risk.
            - 'WR': Worst Realization (Minimax).
            - 'RG': Range of returns.
            - 'CVRG': CVaR range of returns.
            - 'TGRG': Tail Gini range of returns.
            - 'MDD': Maximum Drawdown of uncompounded cumulative returns (Calmar Ratio).
            - 'ADD': Average Drawdown of uncompounded cumulative returns.
            - 'CDaR': Conditional Drawdown at Risk of uncompounded cumulative returns.
            - 'EDaR': Entropic Drawdown at Risk of uncompounded cumulative returns.
            - 'RLDaR': Relativistic Drawdown at Risk of uncompounded cumulative returns.
            - 'UCI': Ulcer Index of uncompounded cumulative returns.

        obj : str can be {'MinRisk', 'Utility', 'Sharpe' or 'MaxRet'}.
            Objective function of the optimization model.
            The default is 'Sharpe'. Possible values are:

            - 'MinRisk': Minimize the selected risk measure.
            - 'Utility': Maximize the Utility function :math:`\mu w - l \phi_{i}(w)`.
            - 'Sharpe': Maximize the risk adjusted return ratio based on the selected risk measure.
            - 'MaxRet': Maximize the expected return of the portfolio.

        kelly : str, optional
            Method used to calculate mean return. Possible values are False for
            arithmetic mean return, "approx" for approximate mean logarithmic
            return using first and second moment and "exact" for mean logarithmic
            return. The default is False.
        rf : float, optional
            Risk free rate, must be in the same period of assets returns.
            The default is 0.
        l : scalar, optional
            Risk aversion factor of the 'Utility' objective function.
            The default is 2.
        hist : bool, optional
            Indicate what kind of returns are used to calculate risk measures
            that depends on scenarios (All except 'MV' risk measure).
            If model = 'BL', True means historical covariance and returns and
            False Black Litterman covariance and historical returns.
            If model = 'FM', True means historical covariance and returns and
            False Risk Factor model for covariance and returns.
            If model = 'BL_FM', True means historical covariance and returns,
            False Black Litterman with Risk Factor model for covariance and
            Risk Factor model for returns, and '2' Risk Factor model for
            covariance and returns. The default is True.

        Returns
        -------
        w : DataFrame
            The weights of optimal portfolio.

        """

        # General model Variables

        mu = None
        sigma = None
        returns = None
        if model == "Classic":
            mu = np.array(self.mu, ndmin=2)
            sigma = np.array(self.cov, ndmin=2)
            returns = np.array(self.returns, ndmin=2)
        elif model == "FM":
            mu = np.array(self.mu_fm, ndmin=2)
            if hist == False:
                sigma = np.array(self.cov_fm, ndmin=2)
                returns = np.array(self.returns_fm, ndmin=2)
            elif hist == True:
                sigma = np.array(self.cov, ndmin=2)
                returns = np.array(self.returns, ndmin=2)
        elif model == "BL":
            mu = np.array(self.mu_bl, ndmin=2)
            if hist == False:
                sigma = np.array(self.cov_bl, ndmin=2)
            elif hist == True:
                sigma = np.array(self.cov, ndmin=2)
            returns = np.array(self.returns, ndmin=2)
        elif model == "BL_FM":
            mu = np.array(self.mu_bl_fm, ndmin=2)
            if hist == False:
                sigma = np.array(self.cov_bl_fm, ndmin=2)
                returns = np.array(self.returns_fm, ndmin=2)
            elif hist == True:
                sigma = np.array(self.cov, ndmin=2)
                returns = np.array(self.returns, ndmin=2)
            elif hist == 2:
                sigma = np.array(self.cov_fm, ndmin=2)
                returns = np.array(self.returns_fm, ndmin=2)

        # General Model Variables

        returns = np.array(returns, ndmin=2)
        T, N = returns.shape
        w = cp.Variable((N, 1))
        k = cp.Variable((1, 1))
        rf0 = rf
        gr = cp.Variable((T, 1))

        # SDP Model Variables

        sdpmodel = False
        if rm in ["KT", "SKT"]:
            sdpmodel = True
        elif (
            self.kurt is not None
            or self.skurt is not None
            or self.network_sdp is not None
        ):
            sdpmodel = True
        elif self.upperkt is not None or self.upperskt is not None:
            sdpmodel = True

        if sdpmodel:
            W = cp.Variable((N, N), symmetric=True)
            M1 = cp.vstack([W, w.T])
            if obj == "Sharpe":
                M2 = cp.vstack([w, k])
            else:
                M2 = cp.vstack([w, np.ones((1, 1))])
            M3 = cp.hstack([M1, M2])
            sdpconstraints = [M3 >> 0]

        # MV Model Variables

        penalty_factor = cp.Constant(0)
        if self.network_sdp is None:
            g = cp.Variable(nonneg=True)
            G = sqrtm(sigma)
            risk1 = g**2
            devconstraints = [cp.SOC(g, G.T @ w)]
        elif self.network_sdp is not None:
            risk1 = cp.trace(sigma @ W)
            devconstraints = []
            if rm != "MV" and self.upperdev is None:
                penalty_factor = self.network_penalty * cp.trace(W)

        # Return Variables

        if model == "Classic":
            if kelly == "exact":
                if obj == "Sharpe":
                    ret = 1 / T * cp.sum(gr) - rf0 * k
                else:
                    ret = 1 / T * cp.sum(cp.log(1 + returns @ w))
            elif kelly == "approx":
                if obj == "Sharpe":
                    ret = mu @ w - 0.5 * cp.quad_over_lin(g, k)
                else:
                    ret = mu @ w - 0.5 * g**2
            elif kelly == False:
                ret = mu @ w
        else:
            ret = mu @ w

        # MAD Model Variables

        madmodel = False
        Y = cp.Variable((T, 1))
        u = np.repeat(mu, T, axis=0)
        a = returns - u
        risk2 = cp.sum(Y) / T
        madconstraints = [a @ w * 1000 >= -Y * 1000, Y * 1000 >= 0]

        # Semi Variance Model Variables

        risk3 = cp.norm(Y, "fro") / cp.sqrt(T - 1)

        # CVaR Model Variables

        VaR = cp.Variable((1, 1))
        alpha = self.alpha
        X = returns @ w
        Z1 = cp.Variable((T, 1))
        risk4 = VaR + 1 / (alpha * T) * cp.sum(Z1)
        cvarconstraints = [Z1 * 1000 >= 0, Z1 * 1000 >= -X * 1000 - VaR * 1000]

        # Worst Realization (Minimax) Model Variables

        M = cp.Variable((1, 1))
        risk5 = M
        wrconstraints = [-X <= M]

        # Lower Partial Moment Variables

        lpmmodel = False
        lpm = cp.Variable((T, 1))
        lpmconstraints = [lpm * 1000 >= 0]

        if obj == "Sharpe":
            lpmconstraints += [lpm * 1000 >= rf0 * k * 1000 - X * 1000]
        else:
            lpmconstraints += [lpm * 1000 >= rf0 * 1000 - X * 1000]

        # First Lower Partial Moment (Omega) Model Variables

        risk6 = cp.sum(lpm) / T

        # Second Lower Partial Moment (Sortino) Model Variables

        risk7 = cp.norm(lpm, "fro") / cp.sqrt(T - 1)

        # Drawdown Model Variables

        drawdown = False

        U = cp.Variable((T + 1, 1))
        ddconstraints = [
            U[1:] * 1000 >= U[:-1] * 1000 - X * 1000,
            U[1:] * 1000 >= 0,
            U[0] * 1000 == 0,
        ]

        # Maximum Drawdown Model Variables

        MDD = cp.Variable((1, 1))
        risk8 = MDD
        mddconstraints = [MDD >= U[1:]]

        # Average Drawdown Model Variables

        risk9 = 1 / T * cp.sum(U[1:])

        # Conditional Drawdown Model Variables

        DaR = cp.Variable((1, 1))
        Zd1 = cp.Variable((T, 1))
        risk10 = DaR + 1 / (alpha * T) * cp.sum(Zd1)
        cdarconstraints = [
            Zd1 * 1000 >= U[1:] * 1000 - DaR * 1000,
            Zd1 * 1000 >= 0,
        ]

        # Ulcer Index Model Variables

        risk11 = cp.norm(U[1:] * 1000, "fro") / np.sqrt(T)

        # Entropic Value at Risk Model Variables

        t1 = cp.Variable((1, 1))
        s1 = cp.Variable((1, 1), nonneg=True)
        ui = cp.Variable((T, 1))
        risk12 = t1 + s1 * np.log(1 / (alpha * T))

        if obj == "Sharpe":
            evarconstraints = [cp.sum(ui) * 1000 <= s1 * 1000]
            evarconstraints += [
                cp.constraints.ExpCone(
                    -X * 1000 - t1 * 1000, np.ones((T, 1)) @ s1 * 1000, ui * 1000
                )
            ]
        else:
            evarconstraints = [cp.sum(ui) <= s1]
            evarconstraints += [
                cp.constraints.ExpCone(-X - t1, np.ones((T, 1)) @ s1, ui)
            ]

        # Entropic Drawdown at Risk Model Variables

        t2 = cp.Variable((1, 1))
        s2 = cp.Variable((1, 1), nonneg=True)
        uj = cp.Variable((T, 1))
        risk13 = t2 + s2 * np.log(1 / (alpha * T))

        if obj == "Sharpe":
            edarconstraints = [cp.sum(uj) * 1000 <= s2 * 1000]
            edarconstraints += [
                cp.constraints.ExpCone(
                    U[1:] * 1000 - t2 * 1000,
                    np.ones((T, 1)) @ s2 * 1000,
                    uj * 1000,
                )
            ]
        else:
            edarconstraints = [cp.sum(uj) <= s2]
            edarconstraints += [
                cp.constraints.ExpCone(U[1:] - t2, np.ones((T, 1)) @ s2, uj)
            ]

        # Gini Mean Difference Model Variables

        owamodel = False
        a1 = cp.Variable((T, 1))
        b1 = cp.Variable((T, 1))
        y = cp.Variable((T, 1))
        risk14 = cp.sum(a1 + b1)

        owaconstraints = [returns @ w == y]
        gmd_w = owa.owa_gmd(T) / 2
        onesvec = np.ones((T, 1))
        gmdconstraints = [y @ gmd_w.T <= onesvec @ a1.T + b1 @ onesvec.T]

        # Tail Gini Model Variables

        a2 = cp.Variable((T, 1))
        b2 = cp.Variable((T, 1))
        risk15 = cp.sum(a2 + b2)
        a_sim = self.a_sim

        tg_w = owa.owa_tg(T, alpha=alpha)
        tgconstraints = [y @ tg_w.T <= onesvec @ a2.T + b2 @ onesvec.T]

        # Range Model Variables

        a3 = cp.Variable((T, 1))
        b3 = cp.Variable((T, 1))
        risk16 = cp.sum(a3 + b3)

        rg_w = owa.owa_rg(T)
        rgconstraints = [y @ rg_w.T <= onesvec @ a3.T + b3 @ onesvec.T]

        # CVaR Range Model Variables

        a4 = cp.Variable((T, 1))
        b4 = cp.Variable((T, 1))
        risk17 = cp.sum(a4 + b4)

        if self.beta is None:
            beta = alpha
        else:
            beta = self.beta

        cvrg_w = owa.owa_cvrg(T, alpha=alpha, beta=beta)
        cvrgconstraints = [y @ cvrg_w.T <= onesvec @ a4.T + b4 @ onesvec.T]

        # Tail Gini Range Model Variables

        a5 = cp.Variable((T, 1))
        b5 = cp.Variable((T, 1))
        risk18 = cp.sum(a5 + b5)

        if self.b_sim is None:
            b_sim = a_sim
        else:
            b_sim = self.b_sim

        tgrg_w = owa.owa_tgrg(T, alpha=alpha, a_sim=a_sim, beta=beta, b_sim=b_sim)
        tgrgconstraints = [y @ tgrg_w.T <= onesvec @ a5.T + b5 @ onesvec.T]

        # Kurtosis Model Variables

        if self.kurt is not None:
            ktconstraints = []

            if N > self.n_max_kurt:
                K = 2 * N
                g2 = cp.Variable((K, 1))
                risk19 = cp.pnorm(g2, p=2)
                A = af.block_vec_pq(self.kurt, N, N)
                s_A, V_A = cf.k_eigh(A, K)
                s_A = np.clip(s_A, 0, np.inf)

                Bi = []
                for i in range(K):
                    B = s_A[i] ** 0.5 * V_A[:, i]
                    B = B.reshape((N, N), order="F").real
                    Bi.append(B)

                for i in range(K):
                    ktconstraints += [g2[i, 0] == cp.trace(Bi[i] @ W)]
            else:
                L_2 = self.L_2
                S_2 = self.S_2
                Sqrt_Sigma_4 = S_2 @ self.kurt @ S_2.T
                Sqrt_Sigma_4 = sqrtm(Sqrt_Sigma_4)
                g2 = cp.Variable(nonneg=True)
                risk19 = g2
                z = L_2 @ cp.reshape(cp.vec(W), (N * N, 1))
                ktconstraints += [cp.SOC(g2, Sqrt_Sigma_4 @ z)]

        # Semi Kurtosis Model Variables

        if self.skurt is not None:
            sktconstraints = []

            if N > self.n_max_kurt:
                K = 2 * N
                sg2 = cp.Variable((K, 1))
                risk20 = cp.pnorm(sg2, p=2)
                SA = af.block_vec_pq(self.skurt, N, N)
                s_SA, V_SA = cf.k_eigh(SA, K)
                s_SA = np.clip(s_SA, 0, np.inf)

                SBi = []
                for i in range(K):
                    SB = s_SA[i] ** 0.5 * V_SA[:, i]
                    SB = SB.reshape((N, N), order="F").real
                    SBi.append(SB)

                for i in range(K):
                    sktconstraints += [sg2[i, 0] == cp.trace(SBi[i] @ W)]

            else:
                Sqrt_SSigma_4 = S_2 @ self.skurt @ S_2.T
                Sqrt_SSigma_4 = sqrtm(Sqrt_SSigma_4)
                sg2 = cp.Variable(nonneg=True)
                risk20 = sg2
                sz = L_2 @ cp.reshape(cp.vec(W), (N * N, 1))
                sktconstraints += [cp.SOC(sg2, Sqrt_SSigma_4 @ sz)]

        # Relativistic Value at Risk Variables

        kappa = self.kappa
        t3 = cp.Variable((1, 1))
        s3 = cp.Variable((1, 1), nonneg=True)
        omega3 = cp.Variable((T, 1))
        psi3 = cp.Variable((T, 1))
        theta3 = cp.Variable((T, 1))
        epsilon3 = cp.Variable((T, 1))

        rlvarconstraints = [
            cp.constraints.power.PowCone3D(
                s3 * (1 + kappa) / (2 * kappa) * onesvec,
                psi3 * (1 + kappa) / kappa,
                epsilon3,
                1 / (1 + kappa),
            ),
            cp.constraints.power.PowCone3D(
                omega3 / (1 - kappa),
                theta3 / kappa,
                -s3 / (2 * kappa) * onesvec,
                (1 - kappa),
            ),
            -X * 1000 - t3 * 1000 + epsilon3 * 1000 + omega3 * 1000 <= 0,
        ]

        ln_k = ((1 / (alpha * T)) ** kappa - (1 / (alpha * T)) ** (-kappa)) / (
            2 * kappa
        )
        risk21 = t3 + ln_k * s3 + cp.sum(psi3 + theta3)

        # Relativistic Drawdown at Risk Variables

        t4 = cp.Variable((1, 1))
        s4 = cp.Variable((1, 1), nonneg=True)
        omega4 = cp.Variable((T, 1))
        psi4 = cp.Variable((T, 1))
        theta4 = cp.Variable((T, 1))
        epsilon4 = cp.Variable((T, 1))

        rldarconstraints = [
            cp.constraints.power.PowCone3D(
                s4 * (1 + kappa) / (2 * kappa) * onesvec,
                psi4 * (1 + kappa) / kappa,
                epsilon4,
                1 / (1 + kappa),
            ),
            cp.constraints.power.PowCone3D(
                omega4 / (1 - kappa),
                theta4 / kappa,
                -s4 / (2 * kappa) * onesvec,
                (1 - kappa),
            ),
            U[1:] * 1000 - t4 * 1000 + epsilon4 * 1000 + omega4 * 1000 <= 0,
        ]

        risk22 = t4 + ln_k * s4 + cp.sum(psi4 + theta4)

        # Cardinality Boolean Variables

        flag_int = False
        if self.card is not None or self.network_ip is not None:
            flag_int = True
            if obj == "Sharpe":
                e = cp.Variable((mu.shape[1], 1), boolean=True)
                e1 = cp.Variable((mu.shape[1], 1))
            else:
                e = cp.Variable((mu.shape[1], 1), boolean=True)

        # Problem Weight Constraints

        if obj == "Sharpe":
            constraints = [cp.sum(w) == self.budget * k, k * 1000 >= 0]
            if self.sht == False:
                constraints += [w <= self.upperlng * k, w * 1000 >= 0]
                if flag_int:
                    constraints += [
                        e1 <= k,
                        e1 >= 0,
                        e1 <= 100000 * e,
                        e1 >= k - 100000 * (1 - e),
                        w <= self.upperlng * e1,
                    ]
                    # Cardinality Constraint
                    if self.card is not None:
                        constraints += [
                            cp.sum(e) <= self.card,
                        ]
                    # Network IP Constraint
                    if self.network_ip is not None:
                        constraints += [
                            np.unique(self.network_ip + np.identity(N), axis=0) @ e <= 1
                        ]

            elif self.sht == True:
                constraints += [
                    cp.sum(cp.pos(w)) * 1000 <= self.upperlng * k * 1000,
                    cp.sum(cp.neg(w)) * 1000 <= self.uppersht * k * 1000,
                ]
                if flag_int:
                    constraints += [
                        e1 <= k,
                        e1 >= 0,
                        e1 <= 100000 * e,
                        e1 >= k - 100000 * (1 - e),
                        w >= -self.uppersht * e1,
                        w <= self.upperlng * e1,
                    ]
                    # Cardinality Constraint
                    if self.card is not None:
                        constraints += [
                            cp.sum(e) <= self.card,
                        ]
                    # Network IP Constraint
                    if self.network_ip is not None:
                        constraints += [
                            np.unique(self.network_ip + np.identity(N), axis=0) @ e <= 1
                        ]

        else:
            constraints = [cp.sum(w) == self.budget]
            if self.sht == False:
                constraints += [w <= self.upperlng, w * 1000 >= 0]
                if flag_int:
                    constraints += [
                        w <= self.upperlng * e,
                    ]
                    # Cardinality Constraint
                    if self.card is not None:
                        constraints += [
                            cp.sum(e) <= self.card,
                        ]
                    # Network IP Constraint
                    if self.network_ip is not None:
                        constraints += [
                            np.unique(self.network_ip + np.identity(N), axis=0) @ e <= 1
                        ]

            elif self.sht == True:
                constraints += [
                    cp.sum(cp.pos(w)) * 1000 <= self.upperlng * 1000,
                    cp.sum(cp.neg(w)) * 1000 <= self.uppersht * 1000,
                ]
                if flag_int:
                    constraints += [
                        w >= -self.uppersht * e,
                        w <= self.upperlng * e,
                    ]
                    # Cardinality Constraint
                    if self.card is not None:
                        constraints += [
                            cp.sum(e) <= self.card,
                        ]
                    # Network IP Constraint
                    if self.network_ip is not None:
                        constraints += [
                            np.unique(self.network_ip + np.identity(N), axis=0) @ e <= 1
                        ]

        # Problem Linear Constraints

        if self.ainequality is not None and self.binequality is not None:
            A = np.array(self.ainequality, ndmin=2) * 1000
            B = np.array(self.binequality, ndmin=2) * 1000
            if obj == "Sharpe":
                constraints += [A @ w - B @ k >= 0]
            else:
                constraints += [A @ w - B >= 0]

        # Number of Effective Assets Constraints

        if self.nea is not None:
            if obj == "Sharpe":
                constraints += [cp.norm(w, "fro") <= 1 / self.nea**0.5 * k]
            else:
                constraints += [cp.norm(w, "fro") <= 1 / self.nea**0.5]

        # Tracking Error Model Variables

        c = np.array(self.benchweights, ndmin=2)
        if self.kindbench == True:
            bench = returns @ c
        elif self.kindbench == False:
            bench = np.array(self.benchindex, ndmin=2)

        # Tracking error Constraints

        if obj == "Sharpe":
            if self.allowTE == True:
                TE_1 = cp.norm(returns @ w - bench @ k, "fro") / cp.sqrt(T - 1)
                constraints += [TE_1 * 1000 <= self.TE * k * 1000]
        else:
            if self.allowTE == True:
                TE_1 = cp.norm(returns @ w - bench, "fro") / cp.sqrt(T - 1)
                constraints += [TE_1 * 1000 <= self.TE * 1000]

        # Turnover Constraints

        if obj == "Sharpe":
            if self.allowTO == True:
                TO_1 = cp.abs(w - c @ k) * 1000
                constraints += [TO_1 <= self.turnover * k * 1000]
        else:
            if self.allowTO == True:
                TO_1 = cp.abs(w - c) * 1000
                constraints += [TO_1 <= self.turnover * 1000]

        # Problem return Constraints

        if self.lowerret is not None:
            if obj == "Sharpe":
                constraints += [ret >= self.lowerret * k]
            else:
                constraints += [ret >= self.lowerret]

        # Problem SDP network constraints

        if self.network_sdp is not None:
            constraints += [cp.multiply(self.network_sdp, W) == 0]

        # Problem centrality measures constraints

        if self.acentrality is not None and self.bcentrality is not None:
            if obj == "Sharpe":
                constraints += [self.acentrality @ w == self.bcentrality * k]
            else:
                constraints += [self.acentrality @ w == self.bcentrality]

        # Problem risk Constraints

        if self.upperdev is not None:
            if self.network_sdp is None:
                if obj == "Sharpe":
                    constraints += [g <= self.upperdev * k]
                else:
                    constraints += [g <= self.upperdev]
                constraints += devconstraints
            else:
                if obj == "Sharpe":
                    constraints += [risk1 <= self.upperdev**2 * k]
                else:
                    constraints += [risk1 <= self.upperdev**2]
                constraints += devconstraints

        if self.uppermad is not None:
            if obj == "Sharpe":
                constraints += [risk2 <= self.uppermad * k / 2]
            else:
                constraints += [risk2 <= self.uppermad / 2]
            madmodel = True

        if self.uppersdev is not None:
            if obj == "Sharpe":
                constraints += [risk3 <= self.uppersdev * k]
            else:
                constraints += [risk3 <= self.uppersdev]
            madmodel = True

        if self.upperCVaR is not None:
            if obj == "Sharpe":
                constraints += [risk4 <= self.upperCVaR * k]
            else:
                constraints += [risk4 <= self.upperCVaR]
            constraints += cvarconstraints

        if self.upperwr is not None:
            if obj == "Sharpe":
                constraints += [-X <= self.upperwr * k]
            else:
                constraints += [-X <= self.upperwr]
            constraints += wrconstraints

        if self.upperflpm is not None:
            if obj == "Sharpe":
                constraints += [risk6 <= self.upperflpm * k]
            else:
                constraints += [risk6 <= self.upperflpm]
            lpmmodel = True

        if self.upperslpm is not None:
            if obj == "Sharpe":
                constraints += [risk7 <= self.upperslpm * k]
            else:
                constraints += [risk7 <= self.upperslpm]
            lpmmodel = True

        if self.uppermdd is not None:
            if obj == "Sharpe":
                constraints += [U[1:] <= self.uppermdd * k]
            else:
                constraints += [U[1:] <= self.uppermdd]
            constraints += mddconstraints
            drawdown = True

        if self.upperadd is not None:
            if obj == "Sharpe":
                constraints += [risk9 <= self.upperadd * k]
            else:
                constraints += [risk9 <= self.upperadd]
            drawdown = True

        if self.upperCDaR is not None:
            if obj == "Sharpe":
                constraints += [risk10 <= self.upperCDaR * k]
            else:
                constraints += [risk10 <= self.upperCDaR]
            constraints += cdarconstraints
            drawdown = True

        if self.upperuci is not None:
            if obj == "Sharpe":
                constraints += [risk11 <= self.upperuci * 1000 * k]
            else:
                constraints += [risk11 <= self.upperuci * 1000]
            drawdown = True

        if self.upperEVaR is not None:
            if obj == "Sharpe":
                constraints += [risk12 <= self.upperEVaR * k]
            else:
                constraints += [risk12 <= self.upperEVaR]
            constraints += evarconstraints

        if self.upperEDaR is not None:
            if obj == "Sharpe":
                constraints += [risk13 <= self.upperEDaR * k]
            else:
                constraints += [risk13 <= self.upperEDaR]
            constraints += edarconstraints

        if self.uppergmd is not None:
            if obj == "Sharpe":
                constraints += [risk14 <= self.uppergmd * k / 2]
            else:
                constraints += [risk14 <= self.uppergmd / 2]
            constraints += gmdconstraints
            owamodel = True

        if self.uppertg is not None:
            if obj == "Sharpe":
                constraints += [risk15 <= self.uppertg * k]
            else:
                constraints += [risk15 <= self.uppertg]
            constraints += tgconstraints
            owamodel = True

        if self.upperrg is not None:
            if obj == "Sharpe":
                constraints += [risk16 <= self.upperrg * k]
            else:
                constraints += [risk16 <= self.upperrg]
            constraints += rgconstraints
            owamodel = True

        if self.uppercvrg is not None:
            if obj == "Sharpe":
                constraints += [risk17 <= self.uppercvrg * k]
            else:
                constraints += [risk17 <= self.uppercvrg]
            constraints += cvrgconstraints
            owamodel = True

        if self.uppertgrg is not None:
            if obj == "Sharpe":
                constraints += [risk18 <= self.uppertgrg * k]
            else:
                constraints += [risk18 <= self.uppertgrg]
            constraints += tgrgconstraints
            owamodel = True

        if self.kurt is not None:
            if self.upperkt is not None:
                if obj == "Sharpe":
                    constraints += [risk19 <= self.upperkt * k]
                else:
                    constraints += [risk19 <= self.upperkt]
                constraints += ktconstraints

        if self.skurt is not None:
            if self.upperskt is not None:
                if obj == "Sharpe":
                    constraints += [risk20 <= self.upperskt * k]
                else:
                    constraints += [risk20 <= self.upperskt]
                constraints += sktconstraints

        if self.upperRLVaR is not None:
            if obj == "Sharpe":
                constraints += [risk21 <= self.upperRLVaR * k]
            else:
                constraints += [risk21 <= self.upperRLVaR]
            constraints += rlvarconstraints

        if self.upperRLDaR is not None:
            if obj == "Sharpe":
                constraints += [risk22 <= self.upperRLDaR * k]
            else:
                constraints += [risk22 <= self.upperRLDaR]
            constraints += rldarconstraints

        # Defining risk function

        if rm == "MV":
            risk = risk1
            if self.upperdev is None:
                constraints += devconstraints
        elif rm == "MAD":
            risk = risk2
            madmodel = True
        elif rm == "MSV":
            risk = risk3
            madmodel = True
        elif rm == "CVaR":
            risk = risk4
            if self.upperCVaR is None:
                constraints += cvarconstraints
        elif rm == "WR":
            risk = risk5
            if self.upperwr is None:
                constraints += wrconstraints
        elif rm == "FLPM":
            risk = risk6
            lpmmodel = True
        elif rm == "SLPM":
            risk = risk7
            lpmmodel = True
        elif rm == "MDD":
            risk = risk8
            drawdown = True
            if self.uppermdd is None:
                constraints += mddconstraints
        elif rm == "ADD":
            risk = risk9
            drawdown = True
        elif rm == "CDaR":
            risk = risk10
            drawdown = True
            if self.upperCDaR is None:
                constraints += cdarconstraints
        elif rm == "UCI":
            risk = risk11
            drawdown = True
            l = l / 1000
        elif rm == "EVaR":
            risk = risk12
            if self.upperEVaR is None:
                constraints += evarconstraints
        elif rm == "EDaR":
            risk = risk13
            drawdown = True
            if self.upperEDaR is None:
                constraints += edarconstraints
        elif rm == "GMD":
            risk = risk14
            owamodel = True
            if self.uppergmd is None:
                constraints += gmdconstraints
        elif rm == "TG":
            risk = risk15
            owamodel = True
            if self.uppertg is None:
                constraints += tgconstraints
        elif rm == "RG":
            risk = risk16
            owamodel = True
            if self.upperrg is None:
                constraints += rgconstraints
        elif rm == "CVRG":
            risk = risk17
            owamodel = True
            if self.uppertgrg is None:
                constraints += cvrgconstraints
        elif rm == "TGRG":
            risk = risk18
            owamodel = True
            if self.uppertg is None:
                constraints += tgrgconstraints
        elif rm == "KT":
            if self.kurt is not None:
                risk = risk19
                if self.upperkt is None:
                    constraints += ktconstraints
            else:
                raise ValueError(
                    "First you need to calculate Cokurtosis Square Matrix."
                )
        elif rm == "SKT":
            if self.skurt is not None:
                risk = risk20
                if self.upperskt is None:
                    constraints += sktconstraints
            else:
                raise ValueError(
                    "First you need to calculate Semi Cokurtosis Square Matrix."
                )
        elif rm == "RLVaR":
            risk = risk21
            if self.upperRLVaR is None:
                constraints += rlvarconstraints
        elif rm == "RLDaR":
            risk = risk22
            drawdown = True
            if self.upperRLDaR is None:
                constraints += rldarconstraints

        if madmodel == True:
            constraints += madconstraints
        if lpmmodel == True:
            constraints += lpmconstraints
        if drawdown == True:
            constraints += ddconstraints
        if owamodel == True:
            constraints += owaconstraints
        if sdpmodel == True:
            constraints += sdpconstraints

        # Frontier Variables

        portafolio = {}

        for i in self.assetslist:
            portafolio.update({i: []})

        # Optimization Process

        # Defining objective function
        if obj == "Sharpe":
            if model == "Classic":
                if kelly == "exact":
                    constraints += [risk <= 1]
                    constraints += [
                        cp.constraints.ExpCone(gr, np.ones((T, 1)) @ k, k + returns @ w)
                    ]
                    objective = cp.Maximize(ret * 1000 - penalty_factor * 1000)
                elif kelly == "approx":
                    constraints += [risk <= 1]
                    if rm != "MV":
                        constraints += devconstraints
                    objective = cp.Maximize(ret - penalty_factor)
                elif kelly == False:
                    constraints += [mu @ w - rf0 * k == 1]
                    objective = cp.Minimize(risk * 1000 + penalty_factor * 1000)
            else:
                constraints += [mu @ w - rf0 * k == 1]
                objective = cp.Minimize(risk * 1000 + penalty_factor * 1000)
        elif obj == "MinRisk":
            objective = cp.Minimize(risk * 1000 + penalty_factor * 1000)
        elif obj == "Utility":
            objective = cp.Maximize(ret - l * risk - penalty_factor)
        elif obj == "MaxRet":
            objective = cp.Maximize(ret * 1000 - penalty_factor * 1000)

        try:
            prob = cp.Problem(objective, constraints)
            for solver in self.solvers:
                try:
                    if len(self.sol_params) == 0:
                        prob.solve(solver=solver)
                    else:
                        prob.solve(solver=solver, **self.sol_params[solver])
                except:
                    pass
                if w.value is not None:
                    break

            if obj == "Sharpe":
                weights = np.array(w.value / k.value, ndmin=2).T
                if rm == "EVaR" or self.upperEVaR is not None:
                    self.z_EVaR = s1.value / k.value
                if rm == "EDaR" or self.upperEDaR is not None:
                    self.z_EDaR = s2.value / k.value
                if rm == "RLVaR" or self.upperRLVaR is not None:
                    self.z_RLVaR = s3.value / k.value
                if rm == "RLDaR" or self.upperRLDaR is not None:
                    self.z_RLDaR = s4.value / k.value
            else:
                weights = np.array(w.value, ndmin=2).T
                if rm == "EVaR" or self.upperEVaR is not None:
                    self.z_EVaR = s1.value
                if rm == "EDaR" or self.upperEDaR is not None:
                    self.z_EDaR = s2.value
                if rm == "RLVaR" or self.upperRLVaR is not None:
                    self.z_RLVaR = s3.value
                if rm == "RLDaR" or self.upperRLDaR is not None:
                    self.z_RLDaR = s4.value

            if self.sht == False:
                weights = np.abs(weights) / np.sum(np.abs(weights)) * self.budget

            for j in self.assetslist:
                portafolio[j].append(weights[0, self.assetslist.index(j)])

        except:
            pass

        try:
            self.optimal = pd.DataFrame(
                portafolio, index=["weights"], dtype=np.float64
            ).T
        except:
            self.optimal = None
            print("The problem doesn't have a solution with actual input parameters")

        return self.optimal

    def rp_optimization(self, model="Classic", rm="MV", rf=0, b=None, hist=True):
        r"""
        This method that calculates the risk parity portfolio using the risk
        budgeting approach :cite:`a-Roncalli` :cite:`a-RichardRoncalli`,
        according to the optimization model selected by the user. The general
        problem that solves is:

        .. math::
            \begin{aligned}
            &\underset{w}{\min} & & \phi(w)\\
            &\text{s.t.} & & b \log(w) \geq c\\
            & & & \mu w \geq \overline{\mu} \\
            & & & Aw \geq B \\
            & & & w \geq 0 \\
            \end{aligned}

        Where:

        :math:`w` are the weights of the portfolio.

        :math:`\mu`: is the vector of expected returns.

        :math:`b` is a vector of risk constraints.

        :math:`Aw \geq B`: is a set of linear constraints.

        :math:`\phi(w)`: is a risk measure.

        :math:`c`: is an arbitrary constant.

        Parameters
        ----------
        model : str can be 'Classic' or 'FM'
            The model used for optimize the portfolio.
            The default is 'Classic'. Possible values are:

            - 'Classic': use estimates of expected return vector and covariance matrix that depends on historical data.
            - 'FM': use estimates of expected return vector and covariance matrix based on a Risk Factor model specified by the user.

        rm : str, optional
            The risk measure used to optimize the portfolio.
            The default is 'MV'. Possible values are:

            - 'MV': Standard Deviation.
            - 'KT': Square Root of Kurtosis.
            - 'MAD': Mean Absolute Deviation.
            - 'GMD': Gini Mean Difference.
            - 'MSV': Semi Standard Deviation.
            - 'SKT': Square Root of Semi Kurtosis.
            - 'FLPM': First Lower Partial Moment (Omega Ratio).
            - 'SLPM': Second Lower Partial Moment (Sortino Ratio).
            - 'CVaR': Conditional Value at Risk.
            - 'TG': Tail Gini.
            - 'EVaR': Entropic Value at Risk.
            - 'RLVaR': Relativistic Value at Risk.
            - 'CVRG': CVaR range of returns.
            - 'TGRG': Tail Gini range of returns.
            - 'CDaR': Conditional Drawdown at Risk of uncompounded cumulative returns.
            - 'EDaR': Entropic Drawdown at Risk of uncompounded cumulative returns.
            - 'RLDaR': Relativistic Drawdown at Risk of uncompounded cumulative returns.
            - 'UCI': Ulcer Index of uncompounded cumulative returns.

        rf : float, optional
            Risk free rate, must be in the same period of assets returns.
            Used for 'FLPM' and 'SLPM'.
            The default is 0.
        b : float, optional
            The vector of risk constraints per asset.
            The default is 1/n (number of assets).
        hist : bool, optional
            Indicate what kind of returns are used to calculate risk measures
            that depends on scenarios (All except 'MV' risk measure).
            If model = 'FM', True means historical covariance and returns and
            False means Risk Factor model for covariance and returns. The
            default is True.

        Returns
        -------
        w : DataFrame
            The weights of optimal portfolio.

        """

        # General model Variables

        mu = None
        sigma = None
        returns = None
        if model == "Classic":
            mu = np.array(self.mu, ndmin=2)
            sigma = np.array(self.cov, ndmin=2)
            returns = np.array(self.returns, ndmin=2)
        elif model == "FM":
            mu = np.array(self.mu_fm, ndmin=2)
            if hist == False:
                sigma = np.array(self.cov_fm, ndmin=2)
                returns = np.array(self.returns_fm, ndmin=2)
            elif hist == True:
                sigma = np.array(self.cov, ndmin=2)
                returns = np.array(self.returns, ndmin=2)

        # General Model Variables

        T, N = returns.shape

        if b is None:
            rb = np.ones((N, 1))
            rb = rb / N
        else:
            self.b = b.copy()
            rb = self.b

        returns = np.array(returns, ndmin=2)
        w = cp.Variable((N, 1))
        k = cp.Variable((1, 1))
        rf0 = rf
        ret = mu @ w
        constraints = []

        # MV Model Variables

        g = cp.Variable(nonneg=True)
        G = sqrtm(sigma)
        risk1 = g**2
        devconstraints = [cp.SOC(g, G.T @ w)]

        # MAD Model Variables

        Y = cp.Variable((T, 1))
        u = np.repeat(mu, T, axis=0)
        a = returns - u
        risk2 = cp.sum(Y) / T
        # madconstraints=[a @ w >= -Y, a @ w <= Y, Y >= 0]
        madconstraints = [a @ w * 1000 >= -Y * 1000, Y * 1000 >= 0]

        # Semi Variance Model Variables

        risk3 = cp.norm(Y, "fro") / cp.sqrt(T - 1)

        # CVaR Model Variables

        VaR = cp.Variable((1, 1))
        alpha = self.alpha
        X = returns @ w
        Z1 = cp.Variable((T, 1))
        risk4 = VaR + 1 / (alpha * T) * cp.sum(Z1)
        cvarconstraints = [Z1 * 1000 >= 0, Z1 * 1000 >= -X * 1000 - VaR * 1000]

        # Lower Partial Moment Variables

        lpm = cp.Variable((T, 1))
        lpmconstraints = [lpm * 1000 >= 0, lpm * 1000 >= rf0 * k * 1000 - X * 1000]

        # First Lower Partial Moment (Omega) Model Variables

        risk6 = cp.sum(lpm) / T

        # Second Lower Partial Moment (Sortino) Model Variables

        risk7 = cp.norm(lpm, "fro") / cp.sqrt(T - 1)

        # Drawdown Model Variables

        # X1 = k + np.cumsum(returns, axis=0) @ w
        U = cp.Variable((T + 1, 1))
        ddconstraints = [
            U[1:] * 1000 >= U[:-1] * 1000 - X * 1000,
            U[1:] * 1000 >= 0,
            U[0] * 1000 == 0,
        ]

        # Conditional Drawdown Model Variables

        DaR = cp.Variable((1, 1))
        Zd = cp.Variable((T, 1))
        risk10 = DaR + 1 / (alpha * T) * cp.sum(Zd)
        cdarconstraints = [
            Zd * 1000 >= U[1:] * 1000 - DaR * 1000,
            Zd * 1000 >= 0,
        ]

        # Ulcer Index Model Variables

        risk11 = cp.norm(U[1:], "fro") / np.sqrt(T)

        # Entropic Value at Risk Model Variables

        t1 = cp.Variable((1, 1))
        s1 = cp.Variable((1, 1), nonneg=True)
        ui = cp.Variable((T, 1))
        risk12 = t1 + s1 * np.log(1 / (alpha * T))
        evarconstraints = [cp.sum(ui) * 1000 <= s1 * 1000]
        evarconstraints += [
            cp.constraints.ExpCone(
                -X * 1000 - t1 * 1000, np.ones((T, 1)) @ s1 * 1000, ui * 1000
            )
        ]

        # Entropic Drawdown at Risk Model Variables

        t2 = cp.Variable((1, 1))
        s2 = cp.Variable((1, 1), nonneg=True)
        uj = cp.Variable((T, 1))
        risk13 = t2 + s2 * np.log(1 / (alpha * T))
        edarconstraints = [cp.sum(uj) * 1000 <= s2 * 1000]
        edarconstraints += [
            cp.constraints.ExpCone(
                U[1:] * 1000 - t2 * 1000,
                np.ones((T, 1)) @ s2 * 1000,
                uj * 1000,
            )
        ]

        # Gini Mean Difference Model Variables

        a1 = cp.Variable((T, 1))
        b1 = cp.Variable((T, 1))
        y = cp.Variable((T, 1))
        risk14 = cp.sum(a1 + b1)

        owaconstraints = [returns @ w == y]
        gmd_w = owa.owa_gmd(T) / 2

        onesvec = np.ones((T, 1))
        gmdconstraints = [
            y @ gmd_w.T * 1000 <= (onesvec @ a1.T + b1 @ onesvec.T) * 1000
        ]

        # Tail Gini Model Variables

        a2 = cp.Variable((T, 1))
        b2 = cp.Variable((T, 1))
        risk15 = cp.sum(a2 + b2)
        a_sim = self.a_sim

        tg_w = owa.owa_tg(T, alpha=alpha, a_sim=a_sim)
        tgconstraints = [y @ tg_w.T * 1000 <= (onesvec @ a2.T + b2 @ onesvec.T) * 1000]

        # CVaR Range Model Variables

        a4 = cp.Variable((T, 1))
        b4 = cp.Variable((T, 1))
        risk17 = cp.sum(a4 + b4)

        if self.beta is None:
            beta = alpha
        else:
            beta = self.beta

        cvrg_w = owa.owa_cvrg(T, alpha=alpha, beta=beta)
        cvrgconstraints = [
            y @ cvrg_w.T * 1000 <= (onesvec @ a4.T + b4 @ onesvec.T) * 1000
        ]

        # Tail Gini Range Model Variables

        a5 = cp.Variable((T, 1))
        b5 = cp.Variable((T, 1))
        risk18 = cp.sum(a5 + b5)

        if self.b_sim is None:
            b_sim = a_sim
        else:
            b_sim = self.b_sim

        tgrg_w = owa.owa_tgrg(T, alpha=alpha, a_sim=a_sim, beta=beta, b_sim=b_sim)
        tgrgconstraints = [
            y @ tgrg_w.T * 1000 <= (onesvec @ a5.T + b5 @ onesvec.T) * 1000
        ]

        # Kurtosis Model Variables

        if self.kurt is not None:
            W = cp.Variable((N, N), symmetric=True)
            M1 = cp.vstack([W, w.T])
            M2 = cp.vstack([w, np.ones((1, 1))])
            M3 = cp.hstack([M1, M2])
            ktconstraints = [M3 >> 0]

            if N > self.n_max_kurt:
                K = 2 * N
                g2 = cp.Variable((K, 1))
                risk19 = cp.pnorm(g2, p=2)
                A = af.block_vec_pq(self.kurt, N, N)
                s_A, V_A = np.linalg.eig(A)
                s_A = np.clip(s_A, 0, np.inf)

                Bi = []
                for i in range(K):
                    B = s_A[i] ** 0.5 * V_A[:, i]
                    B = B.reshape((N, N), order="F").real
                    Bi.append(B)

                for i in range(K):
                    ktconstraints += [g2[i, 0] == cp.trace(Bi[i] @ W)]
            else:
                L_2 = self.L_2
                S_2 = self.S_2
                Sqrt_Sigma_4 = S_2 @ self.kurt @ S_2.T
                Sqrt_Sigma_4 = sqrtm(Sqrt_Sigma_4)
                g2 = cp.Variable(nonneg=True)
                risk19 = g2
                z = L_2 @ cp.reshape(cp.vec(W), (N * N, 1))
                ktconstraints += [cp.SOC(g2, Sqrt_Sigma_4 @ z)]

        # Semi Kurtosis Model Variables

        if self.skurt is not None:
            SW = cp.Variable((N, N), symmetric=True)
            SM1 = cp.vstack([SW, w.T])
            SM2 = cp.vstack([w, np.ones((1, 1))])
            SM3 = cp.hstack([SM1, SM2])
            sktconstraints = [SM3 >> 0]

            if N > self.n_max_kurt:
                K = 2 * N
                sg2 = cp.Variable((K, 1))
                risk20 = cp.pnorm(sg2, p=2)
                SA = af.block_vec_pq(self.skurt, N, N)
                s_SA, V_SA = np.linalg.eig(SA)
                s_SA = np.clip(s_SA, 0, np.inf)

                SBi = []
                for i in range(K):
                    SB = s_SA[i] ** 0.5 * V_SA[:, i]
                    SB = SB.reshape((N, N), order="F").real
                    SBi.append(SB)

                for i in range(K):
                    sktconstraints += [sg2[i, 0] == cp.trace(SBi[i] @ SW)]

            else:
                Sqrt_SSigma_4 = S_2 @ self.skurt @ S_2.T
                Sqrt_SSigma_4 = sqrtm(Sqrt_SSigma_4)
                sg2 = cp.Variable(nonneg=True)
                risk20 = sg2
                sz = L_2 @ cp.reshape(cp.vec(SW), (N * N, 1))
                sktconstraints += [cp.SOC(sg2, Sqrt_SSigma_4 @ sz)]

        # Relativistic Value at Risk Variables

        kappa = self.kappa
        t3 = cp.Variable((1, 1))
        s3 = cp.Variable((1, 1), nonneg=True)
        omega3 = cp.Variable((T, 1))
        psi3 = cp.Variable((T, 1))
        theta3 = cp.Variable((T, 1))
        epsilon3 = cp.Variable((T, 1))

        rlvarconstraints = [
            cp.constraints.power.PowCone3D(
                s3 * (1 + kappa) / (2 * kappa) * onesvec,
                psi3 * (1 + kappa) / kappa,
                epsilon3,
                1 / (1 + kappa),
            ),
            cp.constraints.power.PowCone3D(
                omega3 / (1 - kappa),
                theta3 / kappa,
                -s3 / (2 * kappa) * onesvec,
                (1 - kappa),
            ),
            -X * 1000 - t3 * 1000 + epsilon3 * 1000 + omega3 * 1000 <= 0,
        ]

        ln_k = ((1 / (alpha * T)) ** kappa - (1 / (alpha * T)) ** (-kappa)) / (
            2 * kappa
        )
        risk21 = t3 + ln_k * s3 + onesvec.T @ psi3 + onesvec.T @ theta3

        # Relativistic Drawdown at Risk Variables

        t4 = cp.Variable((1, 1))
        s4 = cp.Variable((1, 1), nonneg=True)
        omega4 = cp.Variable((T, 1))
        psi4 = cp.Variable((T, 1))
        theta4 = cp.Variable((T, 1))
        epsilon4 = cp.Variable((T, 1))

        rldarconstraints = [
            cp.constraints.power.PowCone3D(
                s4 * (1 + kappa) / (2 * kappa) * onesvec,
                psi4 * (1 + kappa) / kappa,
                epsilon4,
                1 / (1 + kappa),
            ),
            cp.constraints.power.PowCone3D(
                omega4 / (1 - kappa),
                theta4 / kappa,
                -s4 / (2 * kappa) * onesvec,
                (1 - kappa),
            ),
            U[1:] * 1000 - t4 * 1000 + epsilon4 * 1000 + omega4 * 1000 <= 0,
        ]

        risk22 = t4 + ln_k * s4 + onesvec.T @ psi4 + onesvec.T @ theta4
        # Problem Linear Constraints

        if self.ainequality is not None and self.binequality is not None:
            A = np.array(self.ainequality, ndmin=2) * 1000
            B = np.array(self.binequality, ndmin=2) * 1000
            constraints += [A @ w - B @ k >= 0]

        # Problem Return Constraint

        if self.lowerret is not None:
            constraints += [ret >= self.lowerret * k]

        # Defining risk function

        if rm == "MV":
            risk = risk1
            constraints += devconstraints
        elif rm == "MAD":
            risk = risk2
            constraints += madconstraints
        elif rm == "MSV":
            risk = risk3
            constraints += madconstraints
        elif rm == "CVaR":
            risk = risk4
            constraints += cvarconstraints
        elif rm == "FLPM":
            risk = risk6
            constraints += lpmconstraints
        elif rm == "SLPM":
            risk = risk7
            constraints += lpmconstraints
        elif rm == "CDaR":
            risk = risk10
            constraints += ddconstraints
            constraints += cdarconstraints
        elif rm == "UCI":
            risk = risk11
            constraints += ddconstraints
        elif rm == "EVaR":
            risk = risk12
            constraints += evarconstraints
        elif rm == "EDaR":
            risk = risk13
            constraints += ddconstraints
            constraints += edarconstraints
        elif rm == "GMD":
            risk = risk14
            constraints += owaconstraints
            constraints += gmdconstraints
        elif rm == "TG":
            risk = risk15
            constraints += owaconstraints
            constraints += tgconstraints
        elif rm == "CVRG":
            risk = risk17
            constraints += owaconstraints
            constraints += cvrgconstraints
        elif rm == "TGRG":
            risk = risk18
            constraints += owaconstraints
            constraints += tgrgconstraints
        elif rm == "KT":
            if self.kurt is not None:
                risk = risk19
                constraints += ktconstraints
            else:
                raise ValueError(
                    "First you need to calculate Cokurtosis Square Matrix."
                )
        elif rm == "SKT":
            if self.skurt is not None:
                risk = risk20
                constraints += sktconstraints
            else:
                raise ValueError(
                    "First you need to calculate Semi Cokurtosis Square Matrix."
                )
        elif rm == "RLVaR":
            risk = risk21
            constraints += rlvarconstraints
        elif rm == "RLDaR":
            risk = risk22
            constraints += ddconstraints
            constraints += rldarconstraints
        # Risk budgeting constraint

        log_w = cp.Variable((N, 1))
        constraints += [
            # rb.T @ cp.log(w) >= 1,
            rb.T @ log_w >= 1,
            cp.ExpCone(log_w * 1000, np.ones((N, 1)) * 1000, w * 1000),
            w * 1000 >= 0,
            cp.sum(w) * 1000 == k * 1000,
        ]

        # Frontier Variables

        portafolio = {}

        for i in self.assetslist:
            portafolio.update({i: []})

        # Optimization Process

        # Defining objective function

        objective = cp.Minimize(risk * 1000)

        try:
            prob = cp.Problem(objective, constraints)
            for solver in self.solvers:
                try:
                    if len(self.sol_params) == 0:
                        prob.solve(solver=solver)
                    else:
                        prob.solve(solver=solver, **self.sol_params[solver])
                except:
                    pass
                if w.value is not None:
                    break

            if rm == "EVaR":
                self.z_EVaR = s1.value
            if rm == "EDaR":
                self.z_EDaR = s2.value
            if rm == "RLVaR":
                self.z_RLVaR = s3.value
            if rm == "RLDaR":
                self.z_RLDaR = s4.value
            weights = np.array(w.value, ndmin=2).T
            weights = np.abs(weights) / np.sum(np.abs(weights))

            for j in self.assetslist:
                portafolio[j].append(weights[0, self.assetslist.index(j)])

        except:
            pass

        try:
            self.rp_optimal = pd.DataFrame(
                portafolio, index=["weights"], dtype=np.float64
            ).T
        except:
            self.rp_optimal = None
            print("The problem doesn't have a solution with actual input parameters")

        return self.rp_optimal

    def rrp_optimization(self, model="Classic", version="A", l=1, b=None, hist=True):
        r"""
        This method that calculates the relaxed risk parity portfolio according
        to the optimization model and version selected by the user
        :cite:`a-GambetaKwon`. The general problem that solves is:

        .. math::
            \begin{aligned}
            &\underset{w}{\min} & & \psi - \gamma & \\
            &\text{s.t.} & & \zeta = \Sigma w \\
            & & & w^{T} \Sigma w \leq \left ( \psi^{2} - \rho^{2} \right ) & \\
            & & & w_{i} \zeta_{i} \geq \gamma^{2} b_{i} & \forall i=1 , \ldots , N \\
            & & & \lambda x^{T} \Theta x \leq \rho^{2} & \\
            & & & \mu w \geq \overline{\mu} & \\
            & & & Aw \geq B & \\
            & & & \sum^{N}_{i=1} w_{i} = 1 & \\
            & & & \psi, \gamma, \rho, w  \geq 0 & \\
            \end{aligned}

        Where:

        :math:`w`: is the vector of weights of the optimum portfolio.

        :math:`\mu`: is the vector of expected returns.

        :math:`\Sigma`: is the covariance matrix of assets returns.

        :math:`\psi`: is the average risk of the portfolio.

        :math:`\gamma`: is the lower bound of each asset risk constribution.

        :math:`b`: is the risk constribution constraint vector.

        :math:`\zeta_{i}`: is the marginal risk of asset :math:`i`.

        :math:`\rho`: is a regularization variable.

        :math:`\lambda`: is a penalty parameter of :math:`\rho`.

        :math:`\Theta = \text{diag}(\Sigma)`

        :math:`Aw \geq B`: is a set of linear constraints.

        Parameters
        ----------
        model : str can be 'Classic' or 'FM'
            The model used for optimize the portfolio.
            The default is 'Classic'. Possible values are:

            - 'Classic': use estimates of expected return vector and covariance matrix that depends on historical data.
            - 'FM': use estimates of expected return vector and covariance matrix based on a Risk Factor model specified by the user.

        version : str can be 'A', 'B' or 'C'
            Relaxed risk parity model version proposed in :cite:`a-RichardRoncalli`.
            The default is 'A'. Possible values are:

            - 'A': without regularization and penalization constraints.
            - 'B': with regularization constraint but without penalization constraint.
            - 'C': with regularization and penalization constraints.

        l : float, optional
            The penalization factor of penalization constraints. Only used with
            version 'C'. The default is 1.
        b : float, optional
            The vector of risk constraints per asset.
            The default is 1/n (number of assets).
        hist : bool, optional
            Indicate what kind of covariance matrix is used.
            If model = 'FM', True means historical covariance and
            False means Risk Factor model for covariance. The default is
            True.

        Returns
        -------
        w : DataFrame
            The weights of optimal portfolio.

        """

        # General model Variables

        mu = None
        sigma = None
        returns = None
        if model == "Classic":
            mu = np.array(self.mu, ndmin=2)
            sigma = np.array(self.cov, ndmin=2)
            returns = np.array(self.returns, ndmin=2)
        elif model == "FM":
            mu = np.array(self.mu_fm, ndmin=2)
            if hist == False:
                sigma = np.array(self.cov_fm, ndmin=2)
                returns = np.array(self.returns_fm, ndmin=2)
            elif hist == True:
                sigma = np.array(self.cov, ndmin=2)
                returns = np.array(self.returns, ndmin=2)

        # General Model Variables

        T, N = returns.shape

        if b is None:
            rb = np.ones((N, 1))
            rb = rb / N
        else:
            self.b = b.copy()
            rb = self.b

        returns = np.array(returns, ndmin=2)
        w = cp.Variable((mu.shape[1], 1))
        ret = mu @ w

        # RRP Model Variables

        G = sqrtm(sigma)
        Theta = np.diag(np.sqrt(np.diag(sigma)))
        psi = cp.Variable(nonneg=True)
        rho = cp.Variable(nonneg=True)
        gamma = cp.Variable(nonneg=True)
        zeta = cp.Variable((mu.shape[1], 1))
        risk = psi - gamma

        # General Model Constraints

        constraints = []
        constraints += [
            zeta == sigma @ w,
            cp.sum(w) == 1,
            gamma >= 0,
            psi >= 0,
            zeta >= 0,
            w >= 0,
        ]

        for i in range(mu.shape[1]):
            constraints += [
                cp.SOC(
                    w[i, 0] + zeta[i, 0],
                    cp.vstack([2 * gamma * rb[i] ** 0.5, w[i, 0] - zeta[i, 0]]),
                )
            ]

        # Specific Model Constraints

        if version == "A":
            constraints += [cp.SOC(psi, G.T @ w)]
        elif version == "B":
            constraints += [
                cp.SOC(
                    2 * psi,
                    cp.vstack([2 * G.T @ w, -2 * rho * np.ones((1, 1))]),
                )
            ]
            constraints += [cp.SOC(rho, G.T @ w)]
        elif version == "C":
            constraints += [
                cp.SOC(
                    2 * psi,
                    cp.vstack([2 * G.T @ w, -2 * rho * np.ones((1, 1))]),
                )
            ]
            constraints += [cp.SOC(rho, l**0.5 * Theta.T @ w)]

        # Problem Linear Constraints

        if self.ainequality is not None and self.binequality is not None:
            A = np.array(self.ainequality, ndmin=2) * 1000
            B = np.array(self.binequality, ndmin=2) * 1000
            constraints += [A @ w - B >= 0]

        # Problem Return Constraint

        if self.lowerret is not None:
            constraints += [ret >= self.lowerret]

        # Frontier Variables

        portafolio = {}

        for i in self.assetslist:
            portafolio.update({i: []})

        # Optimization Process

        # Defining objective function

        objective = cp.Minimize(risk * 1000)

        try:
            prob = cp.Problem(objective, constraints)
            for solver in self.solvers:
                try:
                    if len(self.sol_params) == 0:
                        prob.solve(solver=solver)
                    else:
                        prob.solve(solver=solver, **self.sol_params[solver])
                except:
                    pass
                if w.value is not None:
                    break

            weights = np.array(w.value, ndmin=2).T
            weights = np.abs(weights) / np.sum(np.abs(weights))

            for j in self.assetslist:
                portafolio[j].append(weights[0, self.assetslist.index(j)])

        except:
            pass

        try:
            self.rrp_optimal = pd.DataFrame(
                portafolio, index=["weights"], dtype=np.float64
            ).T
        except:
            self.rrp_optimal = None
            print("The problem doesn't have a solution with actual input parameters")

        return self.rrp_optimal

    def wc_optimization(self, obj="Sharpe", rf=0, l=2, Umu="box", Ucov="box"):
        r"""
        This method that calculates the worst case mean variance portfolio
        according to the objective function and uncertainty sets selected by
        the user.

        Parameters
        ----------
        obj : str can be {'MinRisk', 'Utility', 'Sharpe' or 'MaxRet'}.
            Objective function of the optimization model.
            The default is 'Sharpe'. Possible values are:

            - 'MinRisk': Minimize the worst case formulation of the selected risk measure.
            - 'Utility': Maximize the worst case formulation of the Utility function :math:`\mu w - l \phi_{i}(w)`.
            - 'Sharpe': Maximize the worst case formulation of the risk adjusted return ratio based on the selected risk measure.
            - 'MaxRet': Maximize the worst case formulation of the expected return of the portfolio.

        rf : float, optional
            Risk free rate, must be in the same period of assets returns.
            The default is 0.
        l : scalar, optional
            Risk aversion factor of the 'Utility' objective function.
            The default is 2.
        Umu : str, optional
            The type of uncertainty set for the mean vector used in the model.
            The default is 'box'. Possible values are:

            - 'box': Use a box uncertainty set for the mean vector.
            - 'ellip': Use a elliptical uncertainty set for the mean vector.
            - None: Don't use an uncertainty set for mean vector.

        Ucov : str, optional
            The type of uncertainty set for the covariance matrix used in the model.
            The default is 'box'. Possible values are:

            - 'box': Use a box uncertainty set for the covariance matrix.
            - 'ellip': Use a elliptical uncertainty set for the covariance matrix.
            - None: Don't use an uncertainty set for covariance matrix.

        Returns
        -------
        w : DataFrame
            The weights of optimal portfolio.

        """

        # General model Variables

        mu = self.mu.to_numpy()
        sigma = self.cov.to_numpy()
        returns = self.returns.to_numpy()

        cov_l = self.cov_l.to_numpy()
        cov_u = self.cov_u.to_numpy()
        cov_mu = self.cov_mu.to_numpy()
        cov_sigma = self.cov_sigma.to_numpy()
        d_mu = self.d_mu.to_numpy()
        k_mu = self.k_mu
        k_sigma = self.k_sigma

        T, N = returns.shape
        w = cp.Variable((N, 1))
        Au = cp.Variable((N, N), symmetric=True)
        Al = cp.Variable((N, N), symmetric=True)
        Z = cp.Variable((N, N), symmetric=True)

        k = cp.Variable((1, 1))
        rf0 = rf
        g = cp.Variable(nonneg=True)

        constraints = []

        # Uncertainty Sets for Mean Vector

        if Umu == "box":
            if obj == "Sharpe":
                constraints += [mu @ w - d_mu @ cp.abs(w) - rf0 * k >= 1]
            else:
                ret = mu @ w - d_mu @ cp.abs(w)
        elif Umu == "ellip":
            if obj == "Sharpe":
                constraints += [
                    mu @ w - k_mu * cp.pnorm(sqrtm(cov_mu) @ w, 2) - rf0 * k >= 1
                ]
            else:
                ret = mu @ w - k_mu * cp.pnorm(sqrtm(cov_mu) @ w, 2)
        else:
            if obj == "Sharpe":
                constraints += [mu @ w - rf0 * k >= 1]
            else:
                ret = mu @ w

        # SDP Model Variables

        sdpmodel = False
        if self.network_sdp is not None or Ucov in ["box", "ellip"]:
            sdpmodel = True

        if sdpmodel:
            W = cp.Variable((N, N), symmetric=True)
            M1 = cp.vstack([W, w.T])
            if obj == "Sharpe":
                M2 = cp.vstack([w, k])
            else:
                M2 = cp.vstack([w, np.ones((1, 1))])
            M = cp.hstack([M1, M2])
            sdpconstraints = [M >> 0]

        # Uncertainty Sets for Covariance Matrix

        if Ucov == "box":
            risk = cp.trace(Au @ cov_u) - cp.trace(Al @ cov_l)
            constraints += [W == Au - Al, Au >= 0, Al >= 0]
        elif Ucov == "ellip":
            risk = cp.trace(sigma @ (W + Z))
            risk += k_sigma * cp.pnorm(sqrtm(cov_sigma) @ (cp.vec(W) + cp.vec(Z)), 2)
            constraints += [Z >> 0]
        else:
            if sdpmodel:
                risk = cp.trace(sigma @ w)
            else:
                G = sqrtm(sigma)
                risk = g**2
                constraints += [cp.SOC(g, G.T @ w)]

        # Cardinal Boolean Variables

        flag_int = False
        if self.card is not None or self.network_ip is not None:
            flag_int = True
            if obj == "Sharpe":
                e = cp.Variable((mu.shape[1], 1), boolean=True)
                e1 = cp.Variable((mu.shape[1], 1))
            else:
                e = cp.Variable((mu.shape[1], 1), boolean=True)

        # Problem Weight Constraints

        if obj == "Sharpe":
            constraints += [cp.sum(w) == self.budget * k, k * 1000 >= 0]
            if self.sht == False:
                constraints += [w <= self.upperlng * k, w * 1000 >= 0]
                if flag_int:
                    constraints += [
                        e1 <= k,
                        e1 >= 0,
                        e1 <= 100000 * e,
                        e1 >= k - 100000 * (1 - e),
                        w <= self.upperlng * e1,
                    ]
                    # Cardinality Constraint
                    if self.card is not None:
                        constraints += [
                            cp.sum(e) <= self.card,
                        ]
                    # Network IP Constraint
                    if self.network_ip is not None:
                        constraints += [
                            np.unique(self.network_ip + np.identity(N), axis=0) @ e <= 1
                        ]

            elif self.sht == True:
                constraints += [
                    cp.sum(cp.pos(w)) * 1000 <= self.upperlng * k * 1000,
                    cp.sum(cp.neg(w)) * 1000 <= self.uppersht * k * 1000,
                ]
                if flag_int:
                    constraints += [
                        e1 <= k,
                        e1 >= 0,
                        e1 <= 100000 * e,
                        e1 >= k - 100000 * (1 - e),
                        w >= -self.uppersht * e1,
                        w <= self.upperlng * e1,
                    ]
                    # Cardinality Constraint
                    if self.card is not None:
                        constraints += [
                            cp.sum(e) <= self.card,
                        ]
                    # Network IP Constraint
                    if self.network_ip is not None:
                        constraints += [
                            np.unique(self.network_ip + np.identity(N), axis=0) @ e <= 1
                        ]

        else:
            constraints += [cp.sum(w) == self.budget]
            if self.sht == False:
                constraints += [w <= self.upperlng, w * 1000 >= 0]
                if flag_int:
                    constraints += [
                        w <= self.upperlng * e,
                    ]
                    # Cardinality Constraint
                    if self.card is not None:
                        constraints += [
                            cp.sum(e) <= self.card,
                        ]
                    # Network IP Constraint
                    if self.network_ip is not None:
                        constraints += [
                            np.unique(self.network_ip + np.identity(N), axis=0) @ e <= 1
                        ]

            elif self.sht == True:
                constraints += [
                    cp.sum(cp.pos(w)) * 1000 <= self.upperlng * 1000,
                    cp.sum(cp.neg(w)) * 1000 <= self.uppersht * 1000,
                ]
                if flag_int:
                    constraints += [
                        w >= -self.uppersht * e,
                        w <= self.upperlng * e,
                    ]
                    # Cardinality Constraint
                    if self.card is not None:
                        constraints += [
                            cp.sum(e) <= self.card,
                        ]
                    # Network IP Constraint
                    if self.network_ip is not None:
                        constraints += [
                            np.unique(self.network_ip + np.identity(N), axis=0) @ e <= 1
                        ]

        # Number of effective assets constraints

        if self.nea is not None:
            if obj == "Sharpe":
                constraints += [cp.sum_squares(w) * 1000 <= 1 / self.nea * k * 1000]
            else:
                constraints += [cp.sum_squares(w) * 1000 <= 1 / self.nea * 1000]

        # Tracking Error Model Variables

        c = np.array(self.benchweights, ndmin=2)
        if self.kindbench == True:
            bench = returns @ c
        elif self.kindbench == False:
            bench = np.array(self.benchindex, ndmin=2)

        # Problem additional linear constraints

        if self.ainequality is not None and self.binequality is not None:
            A = np.array(self.ainequality, ndmin=2) * 1000
            B = np.array(self.binequality, ndmin=2) * 1000
            if obj == "Sharpe":
                constraints += [A @ w - B @ k >= 0]
            else:
                constraints += [A @ w - B >= 0]

        # Tracking error Constraints

        if obj == "Sharpe":
            if self.allowTE == True:
                TE_1 = cp.norm(returns @ w - bench @ k, "fro") / np.sqrt(T - 1)
                constraints += [TE_1 <= self.TE * k]
        else:
            if self.allowTE == True:
                TE_1 = cp.norm(returns @ w - bench, "fro") / np.sqrt(T - 1)
                constraints += [TE_1 <= self.TE]

        # Turnover Constraints

        if obj == "Sharpe":
            if self.allowTO == True:
                TO_1 = cp.abs(w - c @ k) * 1000
                constraints += [TO_1 <= self.turnover * k * 1000]
        else:
            if self.allowTO == True:
                TO_1 = cp.abs(w - c) * 1000
                constraints += [TO_1 <= self.turnover * 1000]

        # Problem SDP network constraints

        if self.network_sdp is not None:
            constraints += [cp.multiply(self.network_sdp, W) == 0]

        # Problem centrality measures constraints

        if self.acentrality is not None and self.bcentrality is not None:
            if obj == "Sharpe":
                constraints += [self.acentrality @ w == self.bcentrality * k]
            else:
                constraints += [self.acentrality @ w == self.bcentrality]

        # SDP constraints

        if sdpmodel == True:
            constraints += sdpconstraints

        # Frontier Variables

        portafolio = {}
        for i in self.assetslist:
            portafolio.update({i: []})

        # Optimization Process

        # Defining objective function
        if obj == "Sharpe":
            objective = cp.Minimize(risk * 1000)
        elif obj == "MinRisk":
            objective = cp.Minimize(risk * 1000)
        elif obj == "Utility":
            objective = cp.Maximize(ret * 1000 - l * risk * 1000)
        elif obj == "MaxRet":
            objective = cp.Maximize(ret * 1000)

        try:
            prob = cp.Problem(objective, constraints)
            for solver in self.solvers:
                try:
                    if len(self.sol_params) == 0:
                        prob.solve(solver=solver)
                    else:
                        prob.solve(solver=solver, **self.sol_params[solver])
                except:
                    pass
                if w.value is not None:
                    break

            if obj == "Sharpe":
                weights = np.array(w.value / k.value, ndmin=2).T
            else:
                weights = np.array(w.value, ndmin=2).T

            if self.sht == False:
                weights = np.abs(weights) / np.sum(np.abs(weights)) * self.budget

            for j in self.assetslist:
                portafolio[j].append(weights[0, self.assetslist.index(j)])

        except:
            pass

        try:
            self.wc_optimal = pd.DataFrame(
                portafolio, index=["weights"], dtype=np.float64
            ).T
        except:
            self.wc_optimal = None
            print("The problem doesn't have a solution with actual input parameters")

        return self.wc_optimal

    def owa_optimization(self, obj="Sharpe", owa_w=None, kelly=False, rf=0, l=2):
        r"""
        This method that calculates the owa optimal portfolio according to the
        weight vector given by the user. The general problem that
        solves is:

        .. math::
            \begin{align}
            &\underset{w}{\text{optimize}} & & F(w)\\
            &\text{s. t.} & & Aw \geq B\\
            \end{align}

        Where:

        :math:`F(w)` is the objective function based on an owa risk measure.

        :math:`Aw \geq B` is a set of linear constraints.

        Parameters
        ----------
        obj : str can be {'MinRisk', 'Utility', 'Sharpe' or 'MaxRet'}.
            Objective function of the optimization model.
            The default is 'Sharpe'. Possible values are:

            - 'MinRisk': Minimize the selected risk measure.
            - 'Utility': Maximize the Utility function :math:`\mu w - l \phi_{i}(w)`.
            - 'Sharpe': Maximize the risk adjusted return ratio based on the selected risk measure.

        owa_w : 1darray, optional
            The owa weight used to define the owa risk measure.
            The default is 'MV'. Possible values are:
        kelly : str, optional
            Method used to calculate mean return. Possible values are False for
            arithmetic mean return, "approx" for approximate mean logarithmic
            return using first and second moment and "exact" for mean logarithmic
            return. The default is False.
        rf : float, optional
            Risk free rate, must be in the same period of assets returns.
            The default is 0.
        l : scalar, optional
            Risk aversion factor of the 'Utility' objective function.
            The default is 2.

        Returns
        -------
        w : DataFrame
            The weights of optimal portfolio.

        """

        # General model Variables

        mu = self.mu.to_numpy()
        sigma = self.cov.to_numpy()
        returns = self.returns.to_numpy()
        w = cp.Variable((mu.shape[1], 1))
        k = cp.Variable((1, 1))
        rf0 = rf
        T, N = returns.shape
        gr = cp.Variable((T, 1))

        # MV Model Variables (for approx log returns)

        g = cp.Variable(nonneg=True)
        G = sqrtm(sigma)
        devconstraints = [cp.SOC(g, G.T @ w)]

        # Return Variables

        if kelly == "exact":
            if obj == "Sharpe":
                ret = 1 / T * cp.sum(gr) - rf0 * k
            else:
                ret = 1 / T * cp.sum(cp.log(1 + returns @ w))
        elif kelly == "approx":
            if obj == "Sharpe":
                ret = mu @ w - 0.5 * cp.quad_over_lin(g, k)
            else:
                ret = mu @ w - 0.5 * g**2
        elif kelly == False:
            ret = mu @ w

        # SDP Model Variables

        penalty_factor = cp.Constant(0)
        sdpmodel = False
        if self.network_sdp is not None:
            sdpmodel = True

        if sdpmodel:
            W = cp.Variable((N, N), symmetric=True)
            M1 = cp.vstack([W, w.T])
            if obj == "Sharpe":
                M2 = cp.vstack([w, k])
            else:
                M2 = cp.vstack([w, np.ones((1, 1))])
            M3 = cp.hstack([M1, M2])
            sdpconstraints = [M3 >> 0]
            penalty_factor = self.network_penalty * cp.trace(W)

        # OWA Model Variables

        a = cp.Variable((T, 1))
        b = cp.Variable((T, 1))
        y = cp.Variable((T, 1))
        risk = cp.sum(a + b)

        constraints = []
        constraints += [returns @ w == y]
        if owa_w is None:
            owa_w = owa.owa_gmd(T) / 2

        onesvec = np.ones((T, 1))
        constraints += [y @ owa_w.T <= onesvec @ a.T + b @ onesvec.T]

        # Cardinality Boolean Variables

        flag_int = False
        if self.card is not None or self.network_ip is not None:
            flag_int = True
            if obj == "Sharpe":
                e = cp.Variable((mu.shape[1], 1), boolean=True)
                e1 = cp.Variable((mu.shape[1], 1))
            else:
                e = cp.Variable((mu.shape[1], 1), boolean=True)

        # Problem Weight Constraints

        if obj == "Sharpe":
            constraints += [cp.sum(w) == self.budget * k, k * 1000 >= 0]
            if self.sht == False:
                constraints += [w <= self.upperlng * k, w * 1000 >= 0]
                if flag_int:
                    constraints += [
                        e1 <= k,
                        e1 >= 0,
                        e1 <= 100000 * e,
                        e1 >= k - 100000 * (1 - e),
                        w <= self.upperlng * e1,
                    ]
                    # Cardinality Constraint
                    if self.card is not None:
                        constraints += [
                            cp.sum(e) <= self.card,
                        ]
                    # Network IP Constraint
                    if self.network_ip is not None:
                        constraints += [
                            np.unique(self.network_ip + np.identity(N), axis=0) @ e <= 1
                        ]

            elif self.sht == True:
                constraints += [
                    cp.sum(cp.pos(w)) * 1000 <= self.upperlng * k * 1000,
                    cp.sum(cp.neg(w)) * 1000 <= self.uppersht * k * 1000,
                ]
                if flag_int:
                    constraints += [
                        e1 <= k,
                        e1 >= 0,
                        e1 <= 100000 * e,
                        e1 >= k - 100000 * (1 - e),
                        w >= -self.uppersht * e1,
                        w <= self.upperlng * e1,
                    ]
                    # Cardinality Constraint
                    if self.card is not None:
                        constraints += [
                            cp.sum(e) <= self.card,
                        ]
                    # Network IP Constraint
                    if self.network_ip is not None:
                        constraints += [
                            np.unique(self.network_ip + np.identity(N), axis=0) @ e <= 1
                        ]

        else:
            constraints += [cp.sum(w) == self.budget]
            if self.sht == False:
                constraints += [w <= self.upperlng, w * 1000 >= 0]
                if flag_int:
                    constraints += [
                        w <= self.upperlng * e,
                    ]
                    # Cardinality Constraint
                    if self.card is not None:
                        constraints += [
                            cp.sum(e) <= self.card,
                        ]
                    # Network IP Constraint
                    if self.network_ip is not None:
                        constraints += [
                            np.unique(self.network_ip + np.identity(N), axis=0) @ e <= 1
                        ]

            elif self.sht == True:
                constraints += [
                    cp.sum(cp.pos(w)) * 1000 <= self.upperlng * 1000,
                    cp.sum(cp.neg(w)) * 1000 <= self.uppersht * 1000,
                ]
                if flag_int:
                    constraints += [
                        w >= -self.uppersht * e,
                        w <= self.upperlng * e,
                    ]
                    # Cardinality Constraint
                    if self.card is not None:
                        constraints += [
                            cp.sum(e) <= self.card,
                        ]
                    # Network IP Constraint
                    if self.network_ip is not None:
                        constraints += [
                            np.unique(self.network_ip + np.identity(N), axis=0) @ e <= 1
                        ]

        # Problem Linear Constraints

        if self.ainequality is not None and self.binequality is not None:
            A = np.array(self.ainequality, ndmin=2) * 1000
            B = np.array(self.binequality, ndmin=2) * 1000
            if obj == "Sharpe":
                constraints += [A @ w - B @ k >= 0]
            else:
                constraints += [A @ w - B >= 0]

        # Number of Effective Assets Constraints

        if self.nea is not None:
            if obj == "Sharpe":
                constraints += [cp.norm(w, "fro") <= 1 / self.nea**0.5 * k]
            else:
                constraints += [cp.norm(w, "fro") <= 1 / self.nea**0.5]

        # Tracking Error Model Variables

        c = np.array(self.benchweights, ndmin=2)
        if self.kindbench == True:
            bench = returns @ c
        elif self.kindbench == False:
            bench = np.array(self.benchindex, ndmin=2)

        # Tracking error Constraints

        if obj == "Sharpe":
            if self.allowTE == True:
                TE_1 = cp.norm(returns @ w - bench @ k, "fro") / cp.sqrt(T - 1)
                constraints += [TE_1 * 1000 <= self.TE * k * 1000]
        else:
            if self.allowTE == True:
                TE_1 = cp.norm(returns @ w - bench, "fro") / cp.sqrt(T - 1)
                constraints += [TE_1 * 1000 <= self.TE * 1000]

        # Turnover Constraints

        if obj == "Sharpe":
            if self.allowTO == True:
                TO_1 = cp.abs(w - c @ k) * 1000
                constraints += [TO_1 <= self.turnover * k * 1000]
        else:
            if self.allowTO == True:
                TO_1 = cp.abs(w - c) * 1000
                constraints += [TO_1 <= self.turnover * 1000]

        # Problem return Constraints

        if self.lowerret is not None:
            if obj == "Sharpe":
                constraints += [ret >= self.lowerret * k]
            else:
                constraints += [ret >= self.lowerret]

        # Problem SDP network constraints

        if self.network_sdp is not None:
            constraints += [cp.multiply(self.network_sdp, W) == 0]

        # Problem centrality measures constraints

        if self.acentrality is not None and self.bcentrality is not None:
            if obj == "Sharpe":
                constraints += [self.acentrality @ w == self.bcentrality * k]
            else:
                constraints += [self.acentrality @ w == self.bcentrality]

        # SDP constraints

        if sdpmodel == True:
            constraints += sdpconstraints

        # Frontier Variables

        portafolio = {}

        for i in self.assetslist:
            portafolio.update({i: []})

        # Optimization Process

        # Defining objective function
        if obj == "Sharpe":
            if kelly == "exact":
                constraints += [risk <= 1]
                constraints += [cp.ExpCone(gr, np.ones((T, 1)) @ k, k + returns @ w)]
                objective = cp.Maximize(ret * 1000 - penalty_factor * 1000)
            elif kelly == "approx":
                constraints += [risk <= 1]
                constraints += devconstraints
                objective = cp.Maximize(ret * 1000 - penalty_factor * 1000)
            elif kelly == False:
                constraints += [mu @ w - rf0 * k == 1]
                objective = cp.Minimize(risk * 1000 + penalty_factor * 1000)
        elif obj == "MinRisk":
            objective = cp.Minimize(risk + penalty_factor)
        elif obj == "Utility":
            objective = cp.Maximize(ret - l * risk - penalty_factor)
        elif obj == "MaxRet":
            objective = cp.Maximize(ret * 1000 - penalty_factor * 1000)

        try:
            prob = cp.Problem(objective, constraints)
            for solver in self.solvers:
                try:
                    if len(self.sol_params) == 0:
                        prob.solve(solver=solver)
                    else:
                        prob.solve(solver=solver, **self.sol_params[solver])
                except:
                    pass
                if w.value is not None:
                    break

            if obj == "Sharpe":
                weights = np.array(w.value / k.value, ndmin=2).T
            else:
                weights = np.array(w.value, ndmin=2).T

            if self.sht == False:
                weights = np.abs(weights) / np.sum(np.abs(weights)) * self.budget

            for j in self.assetslist:
                portafolio[j].append(weights[0, self.assetslist.index(j)])

        except:
            pass

        try:
            self.owa_optimal = pd.DataFrame(
                portafolio, index=["weights"], dtype=np.float64
            ).T
        except:
            self.owa_optimal = None
            print("The problem doesn't have a solution with actual input parameters")

        return self.owa_optimal

    def frontier_limits(self, model="Classic", rm="MV", kelly=False, rf=0, hist=True):
        r"""
        Method that calculates the minimum risk and maximum return portfolios
        available with current assets and constraints.

        Parameters
        ----------
        model : str, optional
            Methodology used to estimate input parameters.
            The default is 'Classic'.
        rm : str, optional
            The risk measure used to optimize the portfolio.
            The default is 'MV'. Possible values are:

            - 'MV': Standard Deviation.
            - 'KT': Square Root of Kurtosis.
            - 'MAD': Mean Absolute Deviation.
            - 'GMD': Gini Mean Difference.
            - 'MSV': Semi Standard Deviation.
            - 'SKT': Square Root of Semi Kurtosis.
            - 'FLPM': First Lower Partial Moment (Omega Ratio).
            - 'SLPM': Second Lower Partial Moment (Sortino Ratio).
            - 'CVaR': Conditional Value at Risk.
            - 'TG': Tail Gini.
            - 'EVaR': Entropic Value at Risk.
            - 'RLVaR': Relativistic Value at Risk.
            - 'WR': Worst Realization (Minimax).
            - 'RG': Range of returns.
            - 'CVRG': CVaR range of returns.
            - 'TGRG': Tail Gini range of returns.
            - 'MDD': Maximum Drawdown of uncompounded cumulative returns (Calmar Ratio).
            - 'ADD': Average Drawdown of uncompounded cumulative returns.
            - 'CDaR': Conditional Drawdown at Risk of uncompounded cumulative returns.
            - 'EDaR': Entropic Drawdown at Risk of uncompounded cumulative returns.
            - 'RLDaR': Relativistic Drawdown at Risk of uncompounded cumulative returns.
            - 'UCI': Ulcer Index of uncompounded cumulative returns.

        kelly : str, optional
            Method used to calculate mean return. Possible values are False for
            arithmetic mean return, "approx" for approximate mean logarithmic
            return using first and second moment and "exact" for mean logarithmic
            return. The default is False.
        rf : scalar, optional
            Risk free rate. The default is 0.
        hist : bool, optional
            Indicate what kind of returns are used to calculate risk measures
            that depends on scenarios (All except 'MV' risk measure).
            If model = 'BL', True means historical covariance and returns and
            False Black Litterman covariance and historical returns.
            If model = 'FM', True means historical covariance and returns and
            False Risk Factor model for covariance and returns.
            If model = 'BL_FM', True means historical covariance and returns,
            False Black Litterman with Risk Factor model for covariance and
            Risk Factor model for returns, and '2' Risk Factor model for
            covariance and returns. The default is True.

        Returns
        -------
        limits : DataFrame
            A dataframe that containts the weights of the portfolios.

        Notes
        -----
        This method is preferable (faster) to use instead of efficient_frontier
        method to know the range of expected return and expected risk.

        """

        w_min = self.optimization(
            model=model, rm=rm, obj="MinRisk", kelly=kelly, rf=rf, l=0, hist=hist
        )
        w_max = self.optimization(
            model=model, rm=rm, obj="MaxRet", kelly=kelly, rf=rf, l=0, hist=hist
        )

        if w_min is not None and w_max is not None:
            self.limits = pd.concat([w_min, w_max], axis=1)
            self.limits.columns = ["w_min", "w_max"]
            return self.limits
        else:
            raise NameError("The limits of the frontier can't be found")

    def efficient_frontier(
        self,
        model="Classic",
        rm="MV",
        kelly=False,
        points=20,
        rf=0,
        solver='CLARABEL',
        hist=True,
    ):
        r"""
        Method that calculates several portfolios in the efficient frontier
        of the selected risk measure, available with current assets and
        constraints.

        Parameters
        ----------
        model : str, optional
            Methodology used to estimate input parameters.
            The default is 'Classic'.
        rm : str, optional
            The risk measure used to optimize the portfolio.
            The default is 'MV'. Possible values are:

            - 'MV': Standard Deviation.
            - 'KT': Square Root of Kurtosis.
            - 'MAD': Mean Absolute Deviation.
            - 'GMD': Gini Mean Difference.
            - 'MSV': Semi Standard Deviation.
            - 'SKT': Square Root of Semi Kurtosis.
            - 'FLPM': First Lower Partial Moment (Omega Ratio).
            - 'SLPM': Second Lower Partial Moment (Sortino Ratio).
            - 'CVaR': Conditional Value at Risk.
            - 'TG': Tail Gini.
            - 'EVaR': Entropic Value at Risk.
            - 'RLVaR': Relativistic Value at Risk.
            - 'WR': Worst Realization (Minimax).
            - 'RG': Range of returns.
            - 'CVRG': CVaR range of returns.
            - 'TGRG': Tail Gini range of returns.
            - 'MDD': Maximum Drawdown of uncompounded cumulative returns (Calmar Ratio).
            - 'ADD': Average Drawdown of uncompounded cumulative returns.
            - 'CDaR': Conditional Drawdown at Risk of uncompounded cumulative returns.
            - 'EDaR': Entropic Drawdown at Risk of uncompounded cumulative returns.
            - 'RLDaR': Relativistic Drawdown at Risk of uncompounded cumulative returns.
            - 'UCI': Ulcer Index of uncompounded cumulative returns.

        kelly : str, optional
            Method used to calculate mean return. Possible values are False for
            arithmetic mean return, "approx" for approximate mean logarithmic
            return using first and second moment and "exact" for mean logarithmic
            return. The default is False.
        points : scalar, optional
            Number of point calculated from the efficient frontier.
            The default is 50.
        rf : scalar, optional
            Risk free rate. The default is 0.
        solver: str, optional
            Solver available for CVXPY that supports power cone programming. Used to calculate RLVaR and RLDaR.
            The default value is 'CLARABEL'.
        hist : bool, optional
            Indicate what kind of returns are used to calculate risk measures
            that depends on scenarios (All except 'MV' risk measure).
            If model = 'BL', True means historical covariance and returns and
            False Black Litterman covariance and historical returns.
            If model = 'FM', True means historical covariance and returns and
            False Risk Factor model for covariance and returns.
            If model = 'BL_FM', True means historical covariance and returns,
            False Black Litterman with Risk Factor model for covariance and
            Risk Factor model for returns, and '2' Risk Factor model for
            covariance and returns. The default is True.

        Returns
        -------
        frontier : DataFrame
            A dataframe that containts the weights of the portfolios.

        Notes
        -----
        It's recommendable that don't use this method when there are too many
        assets (more than 100) and you are using a scenario based risk measure
        (all except standard deviation). It's preferable to use frontier_limits
        method (faster) to know the range of expected return and expected risk.
        """

        mu = None
        sigma = None
        returns = None
        if model == "Classic":
            mu = np.array(self.mu, ndmin=2)
            sigma = np.array(self.cov, ndmin=2)
            returns = np.array(self.returns, ndmin=2)
        elif model == "FM":
            mu = np.array(self.mu_fm, ndmin=2)
            if hist == False:
                sigma = np.array(self.cov_fm, ndmin=2)
                returns = np.array(self.returns_fm, ndmin=2)
            elif hist == True:
                sigma = np.array(self.cov, ndmin=2)
                returns = np.array(self.returns, ndmin=2)
        elif model == "BL":
            mu = np.array(self.mu_bl, ndmin=2)
            if hist == False:
                sigma = np.array(self.cov_bl, ndmin=2)
            elif hist == True:
                sigma = np.array(self.cov, ndmin=2)
            returns = np.array(self.returns, ndmin=2)
        elif model == "BL_FM":
            mu = np.array(self.mu_bl_fm, ndmin=2)
            if hist == False:
                sigma = np.array(self.cov_bl_fm, ndmin=2)
                returns = np.array(self.returns_fm, ndmin=2)
            elif hist == True:
                sigma = np.array(self.cov, ndmin=2)
                returns = np.array(self.returns, ndmin=2)
            elif hist == 2:
                sigma = np.array(self.cov_fm, ndmin=2)
                returns = np.array(self.returns_fm, ndmin=2)

        alpha = self.alpha
        a_sim = self.a_sim

        if self.beta is None:
            beta = alpha
        else:
            beta = self.beta

        if self.b_sim is None:
            b_sim = a_sim
        else:
            b_sim = self.b_sim

        kappa = self.kappa

        limits = self.frontier_limits(model=model, rm=rm, kelly=kelly, rf=rf, hist=hist)

        w_min = np.array(limits.iloc[:, 0], ndmin=2).T
        w_max = np.array(limits.iloc[:, 1], ndmin=2).T

        ret_min = (mu @ w_min).item()
        ret_max = (mu @ w_max).item()

        if rm == "MV":
            risk_min = np.sqrt(w_min.T @ sigma @ w_min).item()
            risk_max = np.sqrt(w_max.T @ sigma @ w_max).item()
        elif rm == "KT":
            risk_min = rk.Kurtosis(returns @ w_min)
            risk_max = rk.Kurtosis(returns @ w_max)
        elif rm == "MAD":
            risk_min = rk.MAD(returns @ w_min)
            risk_max = rk.MAD(returns @ w_max)
        elif rm == "MSV":
            risk_min = rk.SemiDeviation(returns @ w_min)
            risk_max = rk.SemiDeviation(returns @ w_max)
        elif rm == "SKT":
            risk_min = rk.SemiKurtosis(returns @ w_min)
            risk_max = rk.SemiKurtosis(returns @ w_max)
        elif rm == "CVaR":
            risk_min = rk.CVaR_Hist(returns @ w_min, alpha)
            risk_max = rk.CVaR_Hist(returns @ w_max, alpha)
        elif rm == "WR":
            risk_min = rk.WR(returns @ w_min)
            risk_max = rk.WR(returns @ w_max)
        elif rm == "FLPM":
            risk_min = rk.LPM(returns @ w_min, rf, 1)
            risk_max = rk.LPM(returns @ w_max, rf, 1)
        elif rm == "SLPM":
            risk_min = rk.LPM(returns @ w_min, rf, 2)
            risk_max = rk.LPM(returns @ w_max, rf, 2)
        elif rm == "MDD":
            risk_min = rk.MDD_Abs(returns @ w_min)
            risk_max = rk.MDD_Abs(returns @ w_max)
        elif rm == "ADD":
            risk_min = rk.ADD_Abs(returns @ w_min)
            risk_max = rk.ADD_Abs(returns @ w_max)
        elif rm == "CDaR":
            risk_min = rk.CDaR_Abs(returns @ w_min, alpha)
            risk_max = rk.CDaR_Abs(returns @ w_max, alpha)
        elif rm == "UCI":
            risk_min = rk.UCI_Abs(returns @ w_min)
            risk_max = rk.UCI_Abs(returns @ w_max)
        elif rm == "EVaR":
            risk_min = rk.EVaR_Hist(returns @ w_min, alpha)[0]
            risk_max = rk.EVaR_Hist(returns @ w_max, alpha)[0]
        elif rm == "EDaR":
            risk_min = rk.EDaR_Abs(returns @ w_min, alpha)[0]
            risk_max = rk.EDaR_Abs(returns @ w_max, alpha)[0]
        elif rm == "GMD":
            risk_min = rk.GMD(returns @ w_min)
            risk_max = rk.GMD(returns @ w_max)
        elif rm == "TG":
            risk_min = rk.TG(returns @ w_min, alpha, a_sim)
            risk_max = rk.TG(returns @ w_max, alpha, a_sim)
        elif rm == "RG":
            risk_min = rk.RG(returns @ w_min)
            risk_max = rk.RG(returns @ w_max)
        elif rm == "CVRG":
            risk_min = rk.CVRG(returns @ w_min, alpha, beta)
            risk_max = rk.CVRG(returns @ w_max, alpha, beta)
        elif rm == "TGRG":
            risk_min = rk.TGRG(returns @ w_min, alpha, a_sim, beta, b_sim)
            risk_max = rk.TGRG(returns @ w_max, alpha, a_sim, beta, b_sim)
        elif rm == "RLVaR":
            risk_min = rk.RLVaR_Hist(returns @ w_min, alpha, kappa, solver)
            risk_max = rk.RLVaR_Hist(returns @ w_max, alpha, kappa, solver)
        elif rm == "RLDaR":
            risk_min = rk.RLDaR_Abs(returns @ w_min, alpha, kappa, solver)
            risk_max = rk.RLDaR_Abs(returns @ w_max, alpha, kappa, solver)

        mus = np.linspace(ret_min, ret_max, points)

        risks = np.linspace(risk_min, risk_max, points)

        risk_lims = [
            "upperdev",
            "upperkt",
            "uppermad",
            "uppergmd",
            "uppersdev",
            "upperskt",
            "upperCVaR",
            "uppertg",
            "upperEVaR",
            "upperRLVaR",
            "upperwr",
            "upperrg",
            "uppercvrg",
            "uppertgrg",
            "upperflpm",
            "upperslpm",
            "uppermdd",
            "upperadd",
            "upperCDaR",
            "upperEDaR",
            "upperRLDaR",
            "upperuci",
        ]

        risk_names = [
            "MV",
            "KT",
            "MAD",
            "GMD",
            "MSV",
            "SKT",
            "CVaR",
            "TG",
            "EVaR",
            "RLVaR",
            "WR",
            "RG",
            "CVRG",
            "TGRG",
            "FLPM",
            "SLPM",
            "MDD",
            "ADD",
            "CDaR",
            "EDaR",
            "RLDaR",
            "UCI",
        ]

        item = risk_names.index(rm)

        frontier = []
        n = 0
        for i in range(len(risks)):
            try:
                if n == 0:
                    w = self.optimization(
                        model=model,
                        rm=rm,
                        obj="MinRisk",
                        kelly=kelly,
                        rf=rf,
                        l=0,
                        hist=hist,
                    )
                else:
                    setattr(self, risk_lims[item], risks[i])
                    w = self.optimization(
                        model=model,
                        rm=rm,
                        obj="MaxRet",
                        kelly=kelly,
                        rf=rf,
                        l=0,
                        hist=hist,
                    )
                if w is not None:
                    n += 1
                    frontier.append(w)
            except:
                pass

        setattr(self, risk_lims[item], None)
        self.frontier = pd.concat(frontier, axis=1)
        self.frontier.columns = list(range(n))

        return self.frontier

    def reset_risk_constraints(self):
        r"""
        Reset all risk constraints.

        """
        cons = [
            "lowerret",
            "upperdev",
            "upperkt",
            "uppermad",
            "uppersdev",
            "upperskt",
            "upperCVaR",
            "upperEVaR",
            "upperRLVaR",
            "upperwr",
            "upperflpm",
            "upperslpm",
            "uppermdd",
            "upperadd",
            "upperCDaR",
            "upperEDaR",
            "upperRLDaR",
            "upperuci",
            "uppergmd",
            "uppertg",
            "upperrg",
            "uppercvrg",
            "uppertgrg",
        ]

        for i in cons:
            setattr(self, i, None)

    def reset_linear_constraints(self):
        r"""
        Reset all linear constraints.

        """

        self.ainequality = None
        self.binequality = None
        self.b = None
        self.network_sdp = None
        self.network_penalty = 0.05
        self.network_ip = None
        self.acentrality = None
        self.bcentrality = None

    def reset_inputs(self):
        r"""
        Reset all inputs parameters of optimization models.

        """

        cons = [
            "mu",
            "cov",
            "kurt",
            "skurt",
            "L_2",
            "S_2",
            "mu_fm",
            "cov_fm",
            "mu_bl",
            "cov_bl",
            "mu_bl_fm",
            "cov_bl_fm",
            "returns_fm",
            "cov_l",
            "cov_u",
            "cov_mu",
            "cov_sigma",
            "d_mu",
            "k_mu",
            "k_sigma",
            "z_EVaR",
            "z_EDaR",
            "z_RLVaR",
            "z_RLDaR",
        ]

        for i in cons:
            setattr(self, i, None)

    def reset_all(self):
        r"""
        Reset portfolio object to defatult values.

        """

        self.sht = False
        self.uppersht = 0.2
        self.upperlng = 1
        self.budget = 1
        self.nea = None
        self.card = None
        self._factors = None
        self.B = None
        self.alpha = 0.05
        self.a_sim = 100
        self.beta = None
        self.b_sim = None
        self.kappa = 0.30
        self.n_max_kurt = 50
        self.kindbench = True
        self.benchindex = None
        self.benchweights = None
        self.allowTO = False
        self.turnover = 0.05
        self.allowTE = False
        self.TE = 0.05

        self.reset_risk_constraints()
        self.reset_linear_constraints()
        self.reset_inputs()
