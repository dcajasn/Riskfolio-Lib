import numpy as np
import pandas as pd
import cvxpy as cv
import scipy.stats as st
from scipy.linalg import sqrtm
import riskfolio.RiskFunctions as rk
import riskfolio.ParamsEstimation as pe
import riskfolio.AuxFunctions as af


class Portfolio(object):
    r"""
    Class that creates a portfolio object with all properties needed to
    calculate optimum portfolios.

    Parameters
    ----------
    returns : DataFrame, optional
        A dataframe that containts the returns of the assets.
        The default is None.
    sht : bool, optional
        Indicate if the portfolio consider short positions (negative weights).
        The default is False.
    uppersht : float, optional
        Indicate the maximum value of the sum of short positions.
        The default is 0.2.
    upperlng : float, optional
        Indicate the maximum value of the sum of long positions (positive
        weights). The default is 1.
    factors : DataFrame, optional
        A dataframe that containts the returns of the factors.
        The default is None.
    alpha : float, optional
        Significance level of CVaR and CDaR. The default is 0.01.
    kindbench : bool, optional
        True if the benchmark is a portfolio with detailed weights and False if
        the benchmark is an index. The default is True.
    allowTO : bool, optional
        Indicate if there is turnover constraints. The default is False.
    turnover : float, optional
        The maximum limit of turnover deviatons. The default is 0.05.
    allowTE : bool, optional
        Indicate if there is tracking error constraints.. The default is False.
    TE : float, optional
        The maximum limit of tracking error deviatons. The default is 0.05.
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
    upperdev : float, optional
        Constraint on max level of standard deviation. The default is None.
    uppermad : float, optional
        Constraint on max level of MAD. The default is None.
    uppersdev : float, optional
        Constraint on max level of semi standard deviation. The default is None.
    upperflpm : float, optional
        Constraint on max level of first lower partial moment.
        The default is None.
    upperslpm : float, optional
        Constraint on max level of second lower partial moment.
        The default is None.
    upperCVaR : float, optional
        Constraint on max level of CVaR. The default is None.
    upperEVaR : float, optional
        Constraint on max level of EVaR. The default is None.
    upperwr : float, optional
        Constraint on max level of worst realization. The default is None.
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
        factors=None,
        alpha=0.05,
        kindbench=True,
        allowTO=False,
        turnover=0.05,
        allowTE=False,
        TE=0.05,
        benchindex=None,
        benchweights=None,
        ainequality=None,
        binequality=None,
        upperdev=None,
        uppermad=None,
        uppersdev=None,
        upperflpm=None,
        upperslpm=None,
        upperCVaR=None,
        upperEVaR=None,
        upperwr=None,
        uppermdd=None,
        upperadd=None,
        upperCDaR=None,
        upperEDaR=None,
        upperuci=None,
    ):

        # Optimization Models Options

        self._returns = returns
        self.sht = sht
        self.uppersht = uppersht
        self.upperlng = upperlng
        self._factors = factors
        self.alpha = alpha
        self.kindbench = kindbench
        self.benchindex = benchindex
        self._benchweights = benchweights
        self._ainequality = ainequality
        self._binequality = binequality
        self.upperdev = upperdev
        self.uppermad = uppermad
        self.uppersdev = uppersdev
        self.upperCVaR = upperCVaR
        self.upperEVaR = upperEVaR
        self.upperwr = upperwr
        self.uppermdd = uppermdd
        self.upperadd = upperadd
        self.upperCDaR = upperCDaR
        self.upperEDaR = upperEDaR
        self.upperflpm = upperflpm
        self.upperslpm = upperslpm
        self.upperuci = upperuci
        self.allowTO = allowTO
        self.turnover = turnover
        self.allowTE = allowTE
        self.TE = TE

        # Inputs of Optimization Models

        self.mu = None
        self.cov = None
        self.mu_fm = None
        self.cov_fm = None
        self.mu_bl = None
        self.cov_bl = None
        self.mu_bl_fm = None
        self.cov_bl_fm = None
        self.returns_fm = None
        self.nav_fm = None
        self.z_EVaR = None

        # Inputs of Worst Case Optimization Models

        self.cov_l = None
        self.cov_u = None
        self.cov_mu = None
        self.cov_sigma = None
        self.d_mu = None
        self.k_mu = None
        self.k_sigma = None

        # Solver params

        self.solvers = [cv.ECOS, cv.SCS, cv.OSQP, cv.CVXOPT]
        self.sol_params = {
            # cv.ECOS: {"max_iters": 500, "abstol": 1e-8},
            # cv.SCS: {"max_iters": 2500, "eps": 1e-5},
            # cv.OSQP: {"max_iter": 10000, "eps_abs": 1e-8},
            # cv.CVXOPT: {"max_iters": 500, "abstol": 1e-8},
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
    def nav(self):
        if self._returns is not None and isinstance(self._returns, pd.DataFrame):
            return self._returns.cumsum()

    @property
    def assetslist(self):
        if self._returns is not None and isinstance(self._returns, pd.DataFrame):
            return self._returns.columns.tolist()

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
    def benchweights(self):
        n = self.numassets
        if self._benchweights is not None:
            if self._benchweights.shape[0] == n and self._benchweights.shape[1] == 1:
                a = self._benchweights
            else:
                raise NameError("Weights must have a size of shape (n_assets,1)")
        else:
            a = np.array(np.ones([n, 1]) / n, ndmin=2)
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
            a = np.array(np.ones([n, 1]) / n, ndmin=2)
        self._benchweights = a

    @property
    def ainequality(self):
        a = self._ainequality
        if a is not None:
            if a.shape[1] == self.numassets:
                a = a
            else:
                raise NameError(
                    "The array ainequality must have the same number of columns that assets' number"
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
                    "The matrix ainequality must have the same number of columns that assets' number"
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

    def assets_stats(self, method_mu="hist", method_cov="hist", **kwargs):
        r"""
        Calculate the inputs that will be used by the optimization method when
        we select the input model='Classic'.

        Parameters
        ----------
        method_mu : string
            Method used to estimate mean vector.
            The default is 'hist'.
        method_cov : string
            Method used to estimate covariance matrix.
            The default is 'hist'.
        **kwargs : dict
            All aditional parameters of mean_vector and covar_matrix functions.

        See Also
        --------
        riskfolio.ParamsEstimation.mean_vector
        riskfolio.ParamsEstimation.covar_matrix

        """

        self.mu = pe.mean_vector(self.returns, method=method_mu, **kwargs)
        self.cov = pe.covar_matrix(self.returns, method=method_cov, **kwargs)
        value = af.is_pos_def(self.cov, threshold=1e-8)
        if value == False:
            print("You must convert self.cov to a positive definite matrix")

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
        **kwargs
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
            Indicates if use equilibrum or historical excess returns.
            The default is True.
        **kwargs : dict
            Other variables related to the mean and covariance estimation.

        See Also
        --------
        riskfolio.ParamsEstimation.black_litterman

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
            **kwargs
        )
        self.mu_bl = mu
        self.cov_bl = cov

        value = af.is_pos_def(self.cov_bl, threshold=1e-8)
        if value == False:
            print("You must convert self.cov_bl to a positive definite matrix")

    def factors_stats(self, method_mu="hist", method_cov="hist", **kwargs):
        r"""
        Calculate the inputs that will be used by the optimization method when
        we select the input model='FM'.

        Parameters
        ----------
        method_mu : string
            Method used to estimate mean vector.
            The default is 'hist'.
        method_cov : string
            Method used to estimate covariance matrix.
            The default is 'hist'.
        **kwargs : dict
            All aditional parameters of risk_factors function.

        See Also
        --------
        riskfolio.ParamsEstimation.forward_regression
        riskfolio.ParamsEstimation.backward_regression
        riskfolio.ParamsEstimation.loadings_matrix
        riskfolio.ParamsEstimation.risk_factors

        """
        X = self.factors
        Y = self.returns
        mu, cov, returns, nav = pe.risk_factors(
            X, Y, method_mu=method_mu, method_cov=method_cov, **kwargs
        )

        self.mu_fm = mu
        self.cov_fm = cov
        self.returns_fm = returns
        self.nav_fm = nav

        value = af.is_pos_def(self.cov_fm, threshold=1e-8)
        if value == False:
            print("You must convert self.cov_fm to a positive definite matrix")

    def wc_stats(
        self,
        box="s",
        ellip="s",
        q=0.05,
        n_sim=3000,
        window=3,
        dmu=0.1,
        dcov=0.1,
        seed=0,
    ):
        r"""
        Calculate the inputs that will be used by the wc_optimization method.

        Parameters
        ----------
        box : string
            The method used to estimate the box uncertainty sets. The default is 's'. Posible values are:

            - 's': stationary bootstrapping method.
            - 'c': circular bootstrapping method.
            - 'm': moving bootstrapping method.
            - 'n': assuming normal returns to calculate confidence levels.
            - 'd': delta method, this method increase and decrease by a percentage.

        ellip : string
            The method used to estimate the elliptical uncertainty sets. The default is 's'. Posible values are:

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
        See Also
        --------
        riskfolio.ParamsEstimation.bootstrapping

        """

        if box not in list("scmdn"):
            raise ValueError("box only can be 's', 'c', 'm', 'd' or 'n'")
        if ellip not in list("scmn"):
            raise ValueError("box only can be 's', 'c', 'm' or 'n'")

        X = self.returns
        cols = X.columns.tolist()
        cols_2 = [i + "-" + j for i in cols for j in cols]
        (n, m) = X.shape
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
            mu_u = mu + st.norm.ppf(1 - q / 2) * np.diag(cov) / n ** 2
            mu_l = mu - st.norm.ppf(1 - q / 2) * np.diag(cov) / n ** 2
            d_mu = (mu_u - mu_l) / 2
            d_mu = pd.DataFrame(d_mu, index=[0], columns=cols)

            # Defining confidence level of covariance matrix assuming normal returns
            rs = np.random.RandomState(seed=seed)
            A = st.wishart.rvs(n, cov / n, size=10000, random_state=rs)
            cov_l = np.percentile(A, q=q, axis=0)
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
            cov_mu = cov / n
            cov_mu = np.diag(np.diag(cov_mu))
            cov_mu = pd.DataFrame(cov_mu, index=cols, columns=cols)
            # Covariance of covariance matrix
            K = af.commutation_matrix(cov)
            I = np.identity(m * m)
            cov_sigma = n * (I + K) @ np.kron(cov_mu, cov_mu)
            cov_sigma = np.diag(np.diag(cov_sigma))
            cov_sigma = pd.DataFrame(cov_sigma, index=cols_2, columns=cols_2)

        k_mu = st.chi2.ppf(1 - q, df=m) ** 0.5
        k_sigma = st.chi2.ppf(1 - q, df=m * m) ** 0.5

        self.cov_l = cov_l
        self.cov_u = cov_u
        self.cov_mu = cov_mu
        self.cov_sigma = cov_sigma
        self.d_mu = d_mu
        self.k_mu = k_mu
        self.k_sigma = k_sigma

    def optimization(
        self, model="Classic", rm="MV", obj="Sharpe", rf=0, l=2, hist=True
    ):
        r"""
        This method that calculates the optimum portfolio according to the
        optimization model selected by the user. The general problem that
        solves is:
        
        .. math::
            \begin{align}
            &\underset{x}{\text{optimize}} & & F(w)\\
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
        model : str can be 'Classic', 'BL' or 'FM'
            The model used for optimize the portfolio.
            The default is 'Classic'. Posible values are:

            - 'Classic': use estimates of expected return vector and covariance matrix that depends on historical data.
            - 'BL': use estimates of expected return vector and covariance matrix based on the Black Litterman model.
            - 'FM': use estimates of expected return vector and covariance matrix based on a Risk Factor model specified by the user.
            
        rm : str, optional
            The risk measure used to optimze the portfolio.
            The default is 'MV'. Posible values are:
            
            - 'MV': Standard Deviation.
            - 'MAD': Mean Absolute Deviation.
            - 'MSV': Semi Standard Deviation.
            - 'FLPM': First Lower Partial Moment (Omega Ratio).
            - 'SLPM': Second Lower Partial Moment (Sortino Ratio).
            - 'CVaR': Conditional Value at Risk.
            - 'EVaR': Entropic Value at Risk.
            - 'WR': Worst Realization (Minimax)
            - 'MDD': Maximum Drawdown of uncompounded cumulative returns (Calmar Ratio).
            - 'ADD': Average Drawdown of uncompounded cumulative returns.
            - 'CDaR': Conditional Drawdown at Risk of uncompounded cumulative returns.
            - 'EDaR': Entropic Drawdown at Risk of uncompounded cumulative returns.
            - 'UCI': Ulcer Index of uncompounded cumulative returns.
            
        obj : str can be {'MinRisk', 'Utility', 'Sharpe' or 'MaxRet'}.
            Objective function of the optimization model.
            The default is 'Sharpe'. Posible values are:

            - 'MinRisk': Minimize the selected risk measure.
            - 'Utility': Maximize the Utility function :math:`\mu w - l \phi_{i}(w)`.
            - 'Sharpe': Maximize the risk adjusted return ratio based on the selected risk measure.
            - 'MaxRet': Maximize the expected return of the portfolio.
                
        rf : float, optional
            Risk free rate, must be in the same period of assets returns.
            The default is 0.
        l : scalar, optional
            Risk aversion factor of the 'Utility' objective function.
            The default is 2.
        hist : bool, optional
            Indicate if uses historical or factor estimation of returns to 
            calculate risk measures that depends on scenarios (All except
            'MV' risk measure). The default is True.

        Returns
        -------
        w : DataFrame
            The weights of optimum portfolio.

        """

        # General model Variables

        mu = None
        sigma = None
        returns = None
        if model == "Classic":
            mu = np.array(self.mu, ndmin=2)
            sigma = np.array(self.cov, ndmin=2)
            returns = np.array(self.returns, ndmin=2)
            nav = np.array(self.nav, ndmin=2)
        elif model == "FM":
            mu = np.array(self.mu_fm, ndmin=2)
            if hist == False:
                sigma = np.array(self.cov_fm, ndmin=2)
                returns = np.array(self.returns_fm, ndmin=2)
                nav = np.array(self.nav_fm, ndmin=2)
            elif hist == True:
                sigma = np.array(self.cov, ndmin=2)
                returns = np.array(self.returns, ndmin=2)
                nav = np.array(self.nav, ndmin=2)
        elif model == "BL":
            mu = np.array(self.mu_bl, ndmin=2)
            if hist == False:
                sigma = np.array(self.cov_bl, ndmin=2)
            elif hist == True:
                sigma = np.array(self.cov, ndmin=2)
            returns = np.array(self.returns, ndmin=2)
            nav = np.array(self.nav, ndmin=2)
        elif model == "BL_FM":
            mu = np.array(self.mu_bl_fm, ndmin=2)
            if hist == False:
                sigma = np.array(self.cov_bl_fm, ndmin=2)
                returns = np.array(self.returns_fm, ndmin=2)
                nav = np.array(self.nav_fm, ndmin=2)
            elif hist == True:
                sigma = np.array(self.cov, ndmin=2)
                returns = np.array(self.returns, ndmin=2)
                nav = np.array(self.nav, ndmin=2)

        # General Model Variables

        returns = np.array(returns, ndmin=2)
        w = cv.Variable((mu.shape[1], 1))
        k = cv.Variable((1, 1))
        rf0 = rf
        n = returns.shape[0]
        ret = mu @ w

        # MV Model Variables

        g = cv.Variable(nonneg=True)
        G = np.linalg.cholesky(sigma)
        risk1 = g ** 2
        devconstraints = [cv.SOC(g, G.T @ w)]

        # MAD Model Variables

        madmodel = False
        Y = cv.Variable((returns.shape[0], 1))
        u = np.ones((returns.shape[0], 1)) * mu
        a = returns - u
        risk2 = cv.sum(Y) / n
        # madconstraints=[a @ w >= -Y, a @ w <= Y, Y >= 0]
        madconstraints = [a @ w >= -Y, Y >= 0]

        # Semi Variance Model Variables

        risk3 = cv.norm(Y, "fro") / cv.sqrt(n - 1)

        # CVaR Model Variables

        VaR = cv.Variable((1, 1))
        alpha = self.alpha
        X = returns @ w
        Z = cv.Variable((returns.shape[0], 1))
        risk4 = VaR + 1 / (alpha * n) * cv.sum(Z)
        cvarconstraints = [Z >= 0, Z >= -X - VaR]

        # Worst Realization (Minimax) Model Variables

        M = cv.Variable((1, 1))
        risk5 = M
        wrconstraints = [-X <= M]

        # Lower Partial Moment Variables

        lpmmodel = False
        lpm = cv.Variable((returns.shape[0], 1))
        lpmconstraints = [lpm >= 0]

        if obj == "Sharpe":
            lpmconstraints += [lpm >= rf0 * k - X]
        else:
            lpmconstraints += [lpm >= rf0 - X]

        # First Lower Partial Moment (Omega) Model Variables

        risk6 = cv.sum(lpm) / n

        # Second Lower Partial Moment (Sortino) Model Variables

        risk7 = cv.norm(lpm, "fro") / cv.sqrt(n - 1)

        # Drawdown Model Variables

        drawdown = False
        if obj == "Sharpe":
            X1 = k + nav @ w
        else:
            X1 = 1 + nav @ w

        U = cv.Variable((nav.shape[0] + 1, 1))
        ddconstraints = [U[1:] * 1000 >= X1 * 1000, U[1:] * 1000 >= U[:-1] * 1000]

        if obj == "Sharpe":
            ddconstraints += [U[1:] * 1000 >= k * 1000, U[0] * 1000 == k * 1000]
        else:
            ddconstraints += [U[1:] * 1000 >= 1 * 1000, U[0] * 1000 == 1 * 1000]

        # Maximum Drawdown Model Variables

        MDD = cv.Variable((1, 1))
        risk8 = MDD
        mddconstraints = [MDD >= U[1:] - X1]

        # Average Drawdown Model Variables

        risk9 = 1 / n * cv.sum(U[1:] - X1)

        # Conditional Drawdown Model Variables

        CDaR = cv.Variable((1, 1))
        Zd = cv.Variable((nav.shape[0], 1))
        risk10 = CDaR + 1 / (alpha * n) * cv.sum(Zd)
        cdarconstraints = [
            Zd * 1000 >= U[1:] * 1000 - X1 * 1000 - CDaR * 1000,
            Zd * 1000 >= 0,
        ]

        # Ulcer Index Model Variables

        risk11 = cv.norm(U[1:] * 1000 - X1 * 1000, "fro") / np.sqrt(n)

        # Entropic Value at Risk Model Variables

        t = cv.Variable((1, 1))
        s = cv.Variable((1, 1), nonneg=True)
        ui = cv.Variable((n, 1))
        risk12 = t + s * np.log(1 / (alpha * n))

        if obj == "Sharpe":
            evarconstraints = [cv.sum(ui) * 1000 <= s * 1000]
            evarconstraints += [
                cv.constraints.ExpCone(
                    -X * 1000 - t * 1000, np.ones((n, 1)) @ s * 1000, ui * 1000
                )
            ]
        else:
            evarconstraints = [cv.sum(ui) <= s]
            evarconstraints += [cv.constraints.ExpCone(-X - t, np.ones((n, 1)) @ s, ui)]

        # Entropic Drawdown at Risk Model Variables

        t1 = cv.Variable((1, 1))
        s1 = cv.Variable((1, 1), nonneg=True)
        uj = cv.Variable((n, 1))
        risk13 = t1 + s1 * np.log(1 / (alpha * n))

        if obj == "Sharpe":
            edarconstraints = [cv.sum(uj) * 1000 <= s1 * 1000]
            edarconstraints += [
                cv.constraints.ExpCone(
                    U[1:] * 1000 - X1 * 1000 - t1 * 1000,
                    np.ones((n, 1)) @ s1 * 1000,
                    uj * 1000,
                )
            ]
        else:
            edarconstraints = [cv.sum(uj) <= s1]
            edarconstraints += [
                cv.constraints.ExpCone(U[1:] - X1 - t1, np.ones((n, 1)) @ s1, uj)
            ]

        # Tracking Error Model Variables

        c = np.array(self.benchweights, ndmin=2)
        if self.kindbench == True:
            bench = returns @ c
        elif self.kindbench == False:
            bench = np.array(self.benchindex, ndmin=2)

        # Problem aditional linear constraints

        if obj == "Sharpe":
            constraints = [
                cv.sum(w) == self.upperlng * k,
                k >= 0,
                mu @ w - rf0 * k == 1,
            ]
            if self.sht == False:
                constraints += [w <= self.upperlng * k, w * 1000 >= 0]
            elif self.sht == True:
                constraints += [
                    w <= self.upperlng * k,
                    w >= -self.uppersht * k,
                    cv.sum(cv.neg(w)) <= self.uppersht * k,
                ]
        else:
            constraints = [cv.sum(w) == self.upperlng]
            if self.sht == False:
                constraints += [w <= self.upperlng, w * 1000 >= 0]
            elif self.sht == True:
                constraints += [
                    w <= self.upperlng,
                    w >= -self.uppersht,
                    cv.sum(cv.neg(w)) <= self.uppersht,
                ]

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
                TE_1 = cv.norm(returns @ w - bench @ k, "fro") / cv.sqrt(n - 1)
                constraints += [TE_1 * 1000 <= self.TE * k * 1000]
        else:
            if self.allowTE == True:
                TE_1 = cv.norm(returns @ w - bench, "fro") / cv.sqrt(n - 1)
                constraints += [TE_1 * 1000 <= self.TE * 1000]

        # Turnover Constraints

        if obj == "Sharpe":
            if self.allowTO == True:
                TO_1 = cv.abs(w - c @ k) * 1000
                constraints += [TO_1 <= self.turnover * k * 1000]
        else:
            if self.allowTO == True:
                TO_1 = cv.abs(w - c) * 1000
                constraints += [TO_1 <= self.turnover * 1000]

        # Problem risk Constraints

        if self.upperdev is not None:
            if obj == "Sharpe":
                constraints += [g <= self.upperdev * k]
            else:
                constraints += [g <= self.upperdev]
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
                constraints += [U[1:] - X1 <= self.uppermdd * k]
            else:
                constraints += [U[1:] - X1 <= self.uppermdd]
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

        if madmodel == True:
            constraints += madconstraints
        if lpmmodel == True:
            constraints += lpmconstraints
        if drawdown == True:
            constraints += ddconstraints

        # Frontier Variables

        portafolio = {}

        for i in self.assetslist:
            portafolio.update({i: []})

        # Optimization Process

        # Defining objective function
        if obj == "Sharpe":
            objective = cv.Minimize(risk * 1000)
        elif obj == "MinRisk":
            objective = cv.Minimize(risk * 1000)
        elif obj == "Utility":
            objective = cv.Maximize(ret - l * risk)
        elif obj == "MaxRet":
            objective = cv.Maximize(ret * 1000)

        try:
            prob = cv.Problem(objective, constraints)
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
                if rm == "EVaR":
                    self.z_EVaR = s.value / k.value
            else:
                weights = np.array(w.value, ndmin=2).T
                if rm == "EVaR":
                    self.z_EVaR = s.value

            if self.sht == False:
                weights = np.abs(weights) / np.sum(np.abs(weights))

            for j in self.assetslist:
                portafolio[j].append(weights[0, self.assetslist.index(j)])

        except:
            pass

        try:
            optimum = pd.DataFrame(portafolio, index=["weights"], dtype=np.float64).T
        except:
            optimum = None
            print("The problem doesn't have a solution with actual input parameters")

        return optimum

    def rp_optimization(self, model="Classic", rm="MV", rf=0, b=None, hist=True):
        r"""
        This method that calculates the risk parity portfolio according to the
        optimization model selected by the user. The general problem that
        solves is:
        
        .. math::
            \begin{align}
            &\underset{w}{\min} & & R(w)\\
            &\text{s.t.} & & b \log(w) \geq c\\
            & & & w \geq 0 \\
            \end{align}
        
        Where:
        
        :math:`w` are the weights of the portfolio.
        
        :math:`R(w)` is the risk measure.
    
        :math:`b` is a vector of risk constraints.
        
        Parameters
        ----------
        model : str can be 'Classic' or 'FM'
            The model used for optimize the portfolio.
            The default is 'Classic'. Posible values are:

            - 'Classic': use estimates of expected return vector and covariance matrix that depends on historical data.
            - 'FM': use estimates of expected return vector and covariance matrix based on a Risk Factor model specified by the user.
            
        rm : str, optional
            The risk measure used to optimze the portfolio.
            The default is 'MV'. Posible values are:
            
            - 'MV': Standard Deviation.
            - 'MAD': Mean Absolute Deviation.
            - 'MSV': Semi Standard Deviation.
            - 'FLPM': First Lower Partial Moment (Omega Ratio).
            - 'SLPM': Second Lower Partial Moment (Sortino Ratio).
            - 'CVaR': Conditional Value at Risk.
            - 'EVaR': Entropic Value at Risk.
            - 'CDaR': Conditional Drawdown at Risk of uncompounded cumulative returns.
            - 'EDaR': Entropic Drawdown at Risk of uncompounded cumulative returns.
            - 'UCI': Ulcer Index of uncompounded cumulative returns.

        rf : float, optional
            Risk free rate, must be in the same period of assets returns.
            Used for 'FLPM' and 'SLPM'.
            The default is 0.                
        b : float, optional
            The vector of risk constraints per asset.
            The default is 1/n (number of assets).
        hist : bool, optional
            Indicate if uses historical or factor estimation of returns to 
            calculate risk measures that depends on scenarios (All except
            'MV' risk measure). The default is True.

        Returns
        -------
        w : DataFrame
            The weights of optimum portfolio.

        """

        # General model Variables

        mu = None
        sigma = None
        returns = None
        if model == "Classic":
            mu = np.array(self.mu, ndmin=2)
            sigma = np.array(self.cov, ndmin=2)
            returns = np.array(self.returns, ndmin=2)
            nav = np.array(self.nav, ndmin=2)
        elif model == "FM":
            mu = np.array(self.mu_fm, ndmin=2)
            if hist == False:
                sigma = np.array(self.cov_fm, ndmin=2)
                returns = np.array(self.returns_fm, ndmin=2)
                nav = np.array(self.nav_fm, ndmin=2)
            elif hist == True:
                sigma = np.array(self.cov, ndmin=2)
                returns = np.array(self.returns, ndmin=2)
                nav = np.array(self.nav, ndmin=2)

        # General Model Variables

        if b is None:
            b = np.ones((1, mu.shape[1]))
            b = b / mu.shape[1]

        returns = np.array(returns, ndmin=2)
        w = cv.Variable((mu.shape[1], 1))
        rf0 = rf
        n = returns.shape[0]

        # MV Model Variables

        g = cv.Variable(nonneg=True)
        G = np.linalg.cholesky(sigma)
        risk1 = g ** 2
        devconstraints = [cv.SOC(g, G.T @ w)]

        # MAD Model Variables

        Y = cv.Variable((returns.shape[0], 1))
        u = np.ones((returns.shape[0], 1)) * mu
        a = returns - u
        risk2 = cv.sum(Y) / n
        # madconstraints=[a @ w >= -Y, a @ w <= Y, Y >= 0]
        madconstraints = [a @ w >= -Y, Y >= 0]

        # Semi Variance Model Variables

        risk3 = cv.norm(Y, "fro") / cv.sqrt(n - 1)

        # CVaR Model Variables

        VaR = cv.Variable((1, 1))
        alpha = self.alpha
        X = returns @ w
        Z = cv.Variable((returns.shape[0], 1))
        risk4 = VaR + 1 / (alpha * n) * cv.sum(Z)
        cvarconstraints = [Z >= 0, Z >= -X - VaR]

        # Lower Partial Moment Variables

        lpm = cv.Variable((returns.shape[0], 1))
        lpmconstraints = [lpm >= 0, lpm >= rf0 - X]

        # First Lower Partial Moment (Omega) Model Variables

        risk6 = cv.sum(lpm) / n

        # Second Lower Partial Moment (Sortino) Model Variables

        risk7 = cv.norm(lpm, "fro") / cv.sqrt(n - 1)

        # Drawdown Model Variables

        X1 = 1 + nav @ w
        U = cv.Variable((nav.shape[0] + 1, 1))
        ddconstraints = [
            U[1:] * 1000 >= X1 * 1000,
            U[1:] * 1000 >= U[:-1] * 1000,
            U[1:] * 1000 >= 1 * 1000,
            U[0] * 1000 == 1 * 1000,
        ]

        # Conditional Drawdown Model Variables

        CDaR = cv.Variable((1, 1))
        Zd = cv.Variable((nav.shape[0], 1))
        risk10 = CDaR + 1 / (alpha * n) * cv.sum(Zd)
        cdarconstraints = [
            Zd * 1000 >= U[1:] * 1000 - X1 * 1000 - CDaR * 1000,
            Zd * 1000 >= 0,
        ]

        # Ulcer Index Model Variables

        risk11 = cv.norm(U[1:] - X1, "fro") / np.sqrt(n)

        # Entropic Value at Risk Model Variables

        t = cv.Variable((1, 1))
        s = cv.Variable((1, 1), nonneg=True)
        ui = cv.Variable((n, 1))
        risk12 = t + s * np.log(1 / (alpha * n))
        evarconstraints = [cv.sum(ui) * 1000 <= s * 1000]
        evarconstraints += [
            cv.constraints.ExpCone(
                -X * 1000 - t * 1000, np.ones((n, 1)) @ s * 1000, ui * 1000
            )
        ]

        # Entropic Drawdown at Risk Model Variables

        t1 = cv.Variable((1, 1))
        s1 = cv.Variable((1, 1), nonneg=True)
        uj = cv.Variable((n, 1))
        risk13 = t1 + s1 * np.log(1 / (alpha * n))
        edarconstraints = [cv.sum(uj) * 1000 <= s1 * 1000]
        edarconstraints += [
            cv.constraints.ExpCone(
                U[1:] * 1000 - X1 * 1000 - t1 * 1000,
                np.ones((n, 1)) @ s1 * 1000,
                uj * 1000,
            )
        ]

        # Defining risk function

        constraints = []

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

        # Frontier Variables

        portafolio = {}

        for i in self.assetslist:
            portafolio.update({i: []})

        # Optimization Process

        # Defining objective function

        objective = cv.Minimize(risk * 1000)

        constraints += [b @ cv.log(w) * 1000 >= 1 * 1000, w * 1000 >= 0]

        try:
            prob = cv.Problem(objective, constraints)
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
            rp_optimum = pd.DataFrame(portafolio, index=["weights"], dtype=np.float64).T
        except:
            rp_optimum = None
            print("The problem doesn't have a solution with actual input parameters")

        return rp_optimum

    def wc_optimization(self, obj="Sharpe", rf=0, l=2, Umu="box", Ucov="box"):
        r"""
        This method that calculates the worst case mean variance portfolio
        according to the objective function and uncertainty sets selected by
        the user.

        Parameters
        ----------
        obj : str can be {'MinRisk', 'Utility', 'Sharpe' or 'MaxRet'}.
            Objective function of the optimization model.
            The default is 'Sharpe'. Posible values are:

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
            The default is 'box'. Posible values are:

            - 'box': Use a box uncertainty set for the mean vector.
            - 'ellip': Use a elliptical uncertainty set for the mean vector.
            - None: Don't use an uncertainty set for mean vector.

        Ucov : str, optional
            The type of uncertainty set for the covariance matrix used in the model.
            The default is 'box'. Posible values are:

            - 'box': Use a box uncertainty set for the covariance matrix.
            - 'ellip': Use a elliptical uncertainty set for the covariance matrix.
            - None: Don't use an uncertainty set for covariance matrix.

        Returns
        -------
        w : DataFrame
            The weights of optimum portfolio.

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

        n = mu.shape[1]
        w = cv.Variable((n, 1))
        Au = cv.Variable((n, n), symmetric=True)
        Al = cv.Variable((n, n), symmetric=True)
        X = cv.Variable((n, n), symmetric=True)
        Z = cv.Variable((n, n), symmetric=True)

        k = cv.Variable((1, 1))
        rf0 = rf
        g = cv.Variable(nonneg=True)

        constraints = []

        if Umu == "box":
            if obj == "Sharpe":
                constraints += [mu @ w - d_mu @ cv.abs(w) - rf0 * k >= 1]
            else:
                ret = mu @ w - d_mu @ cv.abs(w)
        elif Umu == "ellip":
            if obj == "Sharpe":
                constraints += [
                    mu @ w - k_mu * cv.norm(sqrtm(cov_mu) @ w) - rf0 * k >= 1
                ]
            else:
                ret = mu @ w - k_mu * cv.norm(sqrtm(cov_mu) @ w)
        else:
            if obj == "Sharpe":
                constraints += [mu @ w - rf0 * k >= 1]
            else:
                ret = mu @ w

        if Ucov == "box":
            M1 = cv.vstack([Au - Al, w.T])
            if obj == "Sharpe":
                M2 = cv.vstack([w, k])
            else:
                M2 = cv.vstack([w, np.ones((1, 1))])
            M = cv.hstack([M1, M2])
            risk = cv.trace(Au @ cov_u) - cv.trace(Al @ cov_l)
            constraints += [M >> 0, Au >= 0, Al >= 0]
        elif Ucov == "ellip":
            M1 = cv.vstack([X, w.T])
            if obj == "Sharpe":
                M2 = cv.vstack([w, k])
            else:
                M2 = cv.vstack([w, np.ones((1, 1))])
            M = cv.hstack([M1, M2])
            risk = cv.trace(sigma @ (X + Z))
            risk += k_sigma * cv.norm(sqrtm(cov_sigma) @ (cv.vec(X) + cv.vec(Z)))
            constraints += [M >> 0, Z >> 0]
        else:
            G = np.linalg.cholesky(sigma)
            risk = g ** 2
            constraints += [cv.SOC(g, G.T @ w)]

        if obj == "Sharpe":
            constraints += [cv.sum(w) == k, k >= 0]
            if self.sht == False:
                constraints += [w <= k * self.upperlng, w >= 0]
            elif self.sht == True:
                constraints += [w <= k * self.upperlng, w >= -k * self.uppersht]
                constraints += [cv.sum(cv.neg(w)) <= k * self.uppersht]
        else:
            constraints += [cv.sum(w) == 1]
            if self.sht == False:
                constraints += [w <= self.upperlng, w >= 0]
            if self.sht == True:
                constraints += [w <= self.upperlng, w >= -self.uppersht]
                constraints += [cv.sum(cv.neg(w)) <= self.uppersht]

        # Tracking Error Model Variables

        c = np.array(self.benchweights, ndmin=2)
        if self.kindbench == True:
            bench = returns @ c
        elif self.kindbench == False:
            bench = np.array(self.benchindex, ndmin=2)

        # Problem aditional linear constraints

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
                TE_1 = cv.norm(returns @ w - bench @ k, "fro") / cv.sqrt(n - 1)
                constraints += [TE_1 <= self.TE * k]
        else:
            if self.allowTE == True:
                TE_1 = cv.norm(returns @ w - bench, "fro") / cv.sqrt(n - 1)
                constraints += [TE_1 <= self.TE]

        # Turnover Constraints

        if obj == "Sharpe":
            if self.allowTO == True:
                TO_1 = cv.abs(w - c @ k) * 1000
                constraints += [TO_1 <= self.turnover * k * 1000]
        else:
            if self.allowTO == True:
                TO_1 = cv.abs(w - c) * 1000
                constraints += [TO_1 <= self.turnover * 1000]

        # Frontier Variables

        portafolio = {}

        for i in self.assetslist:
            portafolio.update({i: []})

        # Optimization Process

        # Defining objective function
        if obj == "Sharpe":
            objective = cv.Minimize(risk)
        elif obj == "MinRisk":
            objective = cv.Minimize(risk)
        elif obj == "Utility":
            objective = cv.Maximize(ret - l * risk)
        elif obj == "MaxRet":
            objective = cv.Maximize(ret)

        try:
            prob = cv.Problem(objective, constraints)
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
            wc_optimum = pd.DataFrame(portafolio, index=["weights"], dtype=np.float64).T
        except:
            wc_optimum = None
            print("The problem doesn't have a solution with actual input parameters")

        return wc_optimum

    def frontier_limits(self, model="Classic", rm="MV", rf=0, hist=True):
        r"""
        Method that calculates the minimum risk and maximum return portfolios
        available with current assets and constraints.

        Parameters
        ----------
        model : str, optional
            Methodology used to estimate input parameters.
            The default is 'Classic'.
        rm : str, optional
            The risk measure used to optimze the portfolio.
            The default is 'MV'. Posible values are:

            - 'MV': Standard Deviation.
            - 'MAD': Mean Absolute Deviation.
            - 'MSV': Semi Standard Deviation.
            - 'FLPM': First Lower Partial Moment (Omega Ratio).
            - 'SLPM': Second Lower Partial Moment (Sortino Ratio).
            - 'CVaR': Conditional Value at Risk.
            - 'EVaR': Entropic Value at Risk.
            - 'WR': Worst Realization (Minimax)
            - 'MDD': Maximum Drawdown of uncompounded cumulative returns (Calmar Ratio).
            - 'ADD': Average Drawdown of uncompounded cumulative returns.
            - 'CDaR': Conditional Drawdown at Risk of uncompounded cumulative returns.
            - 'EDaR': Entropic Drawdown at Risk of uncompounded cumulative returns.
            - 'UCI': Ulcer Index of uncompounded cumulative returns.

        rf : scalar, optional
            Risk free rate. The default is 0.
        hist : bool, optional
            Indicate if uses historical or factor estimation of returns to
            calculate risk measures that depends on scenarios (All except
            'MV' risk measure). The default is True.

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
            model=model, rm=rm, obj="MinRisk", rf=rf, l=0, hist=hist
        )
        w_max = self.optimization(
            model=model, rm=rm, obj="MaxRet", rf=rf, l=0, hist=hist
        )

        if w_min is not None and w_max is not None:
            limits = pd.concat([w_min, w_max], axis=1)
            limits.columns = ["w_min", "w_max"]
            return limits
        else:
            raise NameError("The limits of the frontier can't be found")

    def efficient_frontier(self, model="Classic", rm="MV", points=20, rf=0, hist=True):
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
            The risk measure used to optimze the portfolio.
            The default is 'MV'. Posible values are:

            - 'MV': Standard Deviation.
            - 'MAD': Mean Absolute Deviation.
            - 'MSV': Semi Standard Deviation.
            - 'FLPM': First Lower Partial Moment (Omega Ratio).
            - 'SLPM': Second Lower Partial Moment (Sortino Ratio).
            - 'CVaR': Conditional Value at Risk.
            - 'EVaR': Entropic Value at Risk.
            - 'WR': Worst Realization (Minimax)
            - 'MDD': Maximum Drawdown of uncompounded cumulative returns (Calmar Ratio).
            - 'ADD': Average Drawdown of uncompounded cumulative returns.
            - 'CDaR': Conditional Drawdown at Risk of uncompounded cumulative returns.
            - 'EDaR': Entropic Drawdown at Risk of uncompounded cumulative returns.
            - 'UCI': Ulcer Index of uncompounded cumulative returns.
        
        points : scalar, optional
            Number of point calculated from the efficient frontier.
            The default is 50.
        rf : scalar, optional
            Risk free rate. The default is 0.
        hist : bool, optional
            Indicate if uses historical or factor estimation of returns to
            calculate risk measures that depends on scenarios (All except
            'MV' risk measure). The default is True.

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
            nav = np.array(self.nav, ndmin=2)
        elif model == "FM":
            mu = np.array(self.mu_fm, ndmin=2)
            if hist == False:
                sigma = np.array(self.cov_fm, ndmin=2)
                returns = np.array(self.returns_fm, ndmin=2)
                nav = np.array(self.nav_fm, ndmin=2)
            elif hist == True:
                sigma = np.array(self.cov, ndmin=2)
                returns = np.array(self.returns, ndmin=2)
                nav = np.array(self.nav, ndmin=2)
        elif model == "BL":
            mu = np.array(self.mu_bl, ndmin=2)
            if hist == False:
                sigma = np.array(self.cov_bl, ndmin=2)
            elif hist == True:
                sigma = np.array(self.cov, ndmin=2)
            returns = np.array(self.returns, ndmin=2)
            nav = np.array(self.nav, ndmin=2)
        elif model == "BL_FM":
            mu = np.array(self.mu_bl_fm_2, ndmin=2)
            if hist == False:
                sigma = np.array(self.cov_bl_fm_2, ndmin=2)
                returns = np.array(self.returns_fm, ndmin=2)
                nav = np.array(self.nav_fm, ndmin=2)
            elif hist == True:
                sigma = np.array(self.cov, ndmin=2)
                returns = np.array(self.returns, ndmin=2)
                nav = np.array(self.nav, ndmin=2)

        alpha = self.alpha

        limits = self.frontier_limits(model=model, rm=rm, rf=rf, hist=hist)

        w_min = np.array(limits.iloc[:, 0], ndmin=2).T
        w_max = np.array(limits.iloc[:, 1], ndmin=2).T

        ret_min = (mu @ w_min).item()
        ret_max = (mu @ w_max).item()

        if rm == "MV":
            risk_min = np.sqrt(w_min.T @ sigma @ w_min).item()
            risk_max = np.sqrt(w_max.T @ sigma @ w_max).item()
        elif rm == "MAD":
            risk_min = rk.MAD(returns @ w_min)
            risk_max = rk.MAD(returns @ w_max)
        elif rm == "MSV":
            risk_min = rk.SemiDeviation(returns @ w_min)
            risk_max = rk.SemiDeviation(returns @ w_max)
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

        mus = np.linspace(ret_min, ret_max, points)

        risks = np.linspace(risk_min, risk_max, points)

        risk_lims = [
            "upperdev",
            "uppermad",
            "uppersdev",
            "upperCVaR",
            "upperEVaR",
            "upperwr",
            "upperflpm",
            "upperslpm",
            "uppermdd",
            "upperadd",
            "upperCDaR",
            "upperEDaR",
            "upperuci",
        ]

        risk_names = [
            "MV",
            "MAD",
            "MSV",
            "CVaR",
            "EVaR",
            "WR",
            "FLPM",
            "SLPM",
            "MDD",
            "ADD",
            "CDaR",
            "EDaR",
            "UCI",
        ]

        item = risk_names.index(rm)

        frontier = []
        n = 0
        for i in range(len(risks)):
            try:
                if n == 0:
                    w = self.optimization(
                        model=model, rm=rm, obj="MinRisk", rf=rf, l=0, hist=hist
                    )
                else:
                    setattr(self, risk_lims[item], risks[i])
                    w = self.optimization(
                        model=model, rm=rm, obj="MaxRet", rf=rf, l=0, hist=hist
                    )
                if w is not None:
                    n += 1
                    frontier.append(w)
            except:
                pass

        setattr(self, risk_lims[item], None)
        frontier = pd.concat(frontier, axis=1)
        frontier.columns = list(range(n))

        return frontier

    def reset_risk_constraints(self):
        r"""
        Reset all risk constraints.

        """
        cons = [
            "upperdev",
            "uppermad",
            "uppersdev",
            "upperCVaR",
            "upperEVaR",
            "upperwr",
            "upperflpm",
            "upperslpm",
            "uppermdd",
            "upperadd",
            "upperCDaR",
            "upperEDaR",
            "upperuci",
        ]

        for i in cons:
            setattr(self, i, None)

    def reset_linear_constraints(self):
        r"""
        Reset all linear constraints.

        """

        self.ainequality = None
        self.binequality = None

    def reset_inputs(self):
        r"""
        Reset all inputs parameters of optimization models.

        """

        cons = [
            "mu",
            "cov",
            "mu_fm",
            "cov_fm",
            "mu_bl",
            "cov_bl",
            "mu_bl_fm",
            "cov_bl_fm",
            "returns_fm",
            "nav_fm",
            "cov_l",
            "cov_u",
            "cov_mu",
            "cov_sigma",
            "d_mu",
            "k_mu",
            "k_sigma",
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
        self._factors = None
        self.alpha = 0.01
        self.kindbench = True
        self.benchindex = None
        self._benchweights = None
        self.allowTO = False
        self.turnover = 0.05
        self.allowTE = False
        self.TE = 0.05

        self.reset_risk_constraints()
        self.reset_linear_constraints()
        self.reset_inputs()
