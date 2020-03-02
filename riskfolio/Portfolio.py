import numpy as np
import pandas as pd
import cvxpy as cv
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
        the tracking error constraints are calculated respect to this index.
        The default is None.
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
    upperwr : float, optional
        Constraint on max level of worst realization. The default is None.
    uppermdd : float, optional
        Constraint on max level of maximum drawdown of uncompounded cumulated
        returns. The default is None.
    upperadd : float, optional
        Constraint on max level of average drawdown of uncompounded cumulated
        returns. The default is None.
    upperCDaR : float, optional
        Constraint on max level of conditional drawdown at risk of
        uncompounded cumulated returns. The default is None.

    """

    def __init__(
        self,
        returns=None,
        sht=False,
        uppersht=0.2,
        upperlng=1,
        factors=None,
        alpha=0.01,
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
        upperwr=None,
        uppermdd=None,
        upperadd=None,
        upperCDaR=None,
    ):

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
        self.upperwr = upperwr
        self.uppermdd = uppermdd
        self.upperadd = uppermdd
        self.upperCDaR = upperCDaR
        self.upperflpm = upperflpm
        self.upperslpm = upperslpm
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

    @property
    def returns(self):
        a = self._returns
        if isinstance(a, pd.DataFrame):
            self.assetslist = a.columns.tolist()
            self.nav = a.cumsum()
            self.numassets = len(self.assetslist)
        else:
            raise NameError("returns must be a DataFrame")
        return a

    @returns.setter
    def returns(self, value):
        a = value
        if a is not None and isinstance(a, pd.DataFrame):
            self._returns = a
            self.assetslist = a.columns.tolist()
            self.nav = a.cumsum()
            self.numassets = len(self.assetslist)
        else:
            raise NameError("returns must be a DataFrame")

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
            a = np.matrix(np.ones([n, 1]) / n)
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
            a = np.matrix(np.ones([n, 1]) / n)
        self._benchweights = a

    @property
    def ainequality(self):
        a = self._ainequality
        if a is not None:
            if a.shape[1] == self.numassets:
                a = a
            else:
                raise NameError(
                    "The matrix ainequality must have the same number of columns that assets' number"
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
        Calculate the inputs that will be use by the optimization method when 
        we select the input model='Classic'.

        Parameters
        ----------
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
        Calculate the inputs that will be use by the optimization method when 
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
            w = self.benchweights

        if delta is None:
            a = np.matrix(self.mu) * np.matrix(w)
            delta = (a - rf) / (np.matrix(w).T * np.matrix(self.cov) * np.matrix(w))
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
        Calculate the inputs that will be use by the optimization method when 
        we select the input model='FM'.
        
        Parameters
        ----------
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
            - 'WR': Worst Realization (Minimax)
            - 'MDD': Maximum Drawdown of uncompounded returns (Calmar Ratio).
            - 'ADD': Average Drawdown of uncompounded returns.
            - 'CDaR': Conditional Drawdown at Risk of uncompounded returns.
            
        obj : str can be {'MinRisk', 'Utility', 'Sharpe' or 'MaxRet'.
            Objective function of the optimization model.
            The default is 'Sharpe'. Posible values are:

            - 'MinRisk': Minimize the selected risk measure.
            - 'Utility': Maximize the Utility function :math:`mu w - l \phi_{i}(w)`.
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
            mu = np.matrix(self.mu)
            sigma = np.matrix(self.cov)
            returns = np.matrix(self.returns)
            nav = np.matrix(self.nav)
        elif model == "FM":
            mu = np.matrix(self.mu_fm)
            if hist == False:
                sigma = np.matrix(self.cov_fm)
                returns = np.matrix(self.returns_fm)
                nav = np.matrix(self.nav_fm)
            elif hist == True:
                sigma = np.matrix(self.cov)
                returns = np.matrix(self.returns)
                nav = np.matrix(self.nav)
        elif model == "BL":
            mu = np.matrix(self.mu_bl)
            if hist == False:
                sigma = np.matrix(self.cov_bl)
            elif hist == True:
                sigma = np.matrix(self.cov)
            returns = np.matrix(self.returns)
            nav = np.matrix(self.nav)
        elif model == "BL_FM":
            mu = np.matrix(self.mu_bl_fm)
            if hist == False:
                sigma = np.matrix(self.cov_bl_fm)
                returns = np.matrix(self.returns_fm)
                nav = np.matrix(self.nav_fm)
            elif hist == True:
                sigma = np.matrix(self.cov)
                returns = np.matrix(self.returns)
                nav = np.matrix(self.nav)

        # General Model Variables

        returns = np.matrix(returns)
        w = cv.Variable((mu.shape[1], 1))
        k = cv.Variable((1, 1))
        rf0 = cv.Parameter(nonneg=True)
        rf0.value = rf
        n = cv.Parameter(nonneg=True)
        n.value = returns.shape[0]
        ret = mu * w

        # MV Model Variables

        risk1 = cv.quad_form(w, sigma)
        returns_1 = af.cov_returns(sigma) * 1000
        n1 = cv.Parameter(nonneg=True)
        n1.value = returns_1.shape[0]
        risk1_1 = cv.norm(returns_1 * w, "fro") / cv.sqrt(n1 - 1)

        # MAD Model Variables

        madmodel = False
        Y = cv.Variable((returns.shape[0], 1))
        u = np.matrix(np.ones((returns.shape[0], 1)) * mu)
        a = returns - u
        risk2 = cv.sum(Y) / n
        # madconstraints=[a*w >= -Y, a*w <= Y, Y >= 0]
        madconstraints = [a * w <= Y, Y >= 0]

        # Semi Variance Model Variables

        risk3 = cv.norm(Y, "fro") / cv.sqrt(n - 1)

        # CVaR Model Variables

        alpha1 = self.alpha
        VaR = cv.Variable(1)
        alpha = cv.Parameter(nonneg=True)
        alpha.value = alpha1
        X = returns * w
        Z = cv.Variable((returns.shape[0], 1))
        risk4 = VaR + 1 / (alpha * n) * cv.sum(Z)
        cvarconstraints = [Z >= 0, Z >= -X - VaR]

        # Worst Realization (Minimax) Model Variables

        M = cv.Variable(1)
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
            X1 = k + nav * w
        else:
            X1 = 1 + nav * w

        U = cv.Variable((nav.shape[0] + 1, 1))
        ddconstraints = [U[1:] >= X1, U[1:] >= U[:-1]]

        if obj == "Sharpe":
            ddconstraints += [U[1:] >= k, U[0] == k]
        else:
            ddconstraints += [U[1:] >= 1, U[0] == 1]

        # Maximum Drawdown Model Variables

        MDD = cv.Variable(1)
        risk8 = MDD
        mddconstraints = [MDD >= U[1:] - X1]

        # Average Drawdown Model Variables

        risk9 = 1 / n * cv.sum(U[1:] - X1)

        # Conditional Drawdown Model Variables

        CDaR = cv.Variable(1)
        Zd = cv.Variable((nav.shape[0], 1))
        risk10 = CDaR + 1 / (alpha * n) * cv.sum(Zd)
        cdarconstraints = [Zd >= U[1:] - X1 - CDaR, Zd >= 0]

        # Tracking Error Model Variables

        c = self.benchweights
        if self.kindbench == True:
            bench = np.matrix(returns) * c
        else:
            bench = self.benchindex

        if obj == "Sharpe":
            TE = cv.norm(returns * w - bench * k, "fro") / cv.sqrt(n - 1)
        else:
            TE = cv.norm(returns * w - bench, "fro") / cv.sqrt(n - 1)

        # Problem aditional linear constraints

        if obj == "Sharpe":
            constraints = [
                cv.sum(w) == self.upperlng * k,
                k >= 0,
                mu * w - rf0 * k == 1,
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
            A = np.matrix(self.ainequality)
            B = np.matrix(self.binequality)
            if obj == "Sharpe":
                constraints += [A * w - B * k >= 0]
            else:
                constraints += [A * w - B >= 0]

        # Turnover Constraints

        if obj == "Sharpe":
            if self.allowTO == True:
                constraints += [cv.abs(w - c * k) * 1000 <= self.turnover * k * 1000]
        else:
            if self.allowTO == True:
                constraints += [cv.abs(w - c) * 1000 <= self.turnover * 1000]

        # Tracking error Constraints

        if obj == "Sharpe":
            if self.allowTE == True:
                constraints += [TE <= self.TE * k]
        else:
            if self.allowTE == True:
                constraints += [TE <= self.TE]

        # Problem risk Constraints

        if self.upperdev is not None:
            if obj == "Sharpe":
                constraints += [risk1_1 <= self.upperdev * k]
            else:
                constraints += [risk1 <= self.upperdev ** 2]

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

        # Defining risk function

        if rm == "MV":
            if model != "Classic":
                risk = risk1_1
            elif model == "Classic":
                risk = risk1
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

        # Defining solvers
        solvers = [cv.ECOS, cv.SCS, cv.OSQP, cv.CVXOPT, cv.GLPK]

        # Defining objective function
        if obj == "Sharpe":
            if rm != "Classic":
                objective = cv.Minimize(risk)
            elif rm == "Classic":
                objective = cv.Minimize(risk * 1000)
        elif obj == "MinRisk":
            objective = cv.Minimize(risk)
        elif obj == "Utility":
            objective = cv.Maximize(ret - l * risk)
        elif obj == "MaxRet":
            objective = cv.Maximize(ret)

        try:
            prob = cv.Problem(objective, constraints)
            for solver in solvers:
                try:
                    prob.solve(
                        solver=solver, parallel=True, max_iters=2000, abstol=1e-10
                    )
                except:
                    pass
                if w.value is not None:
                    break

            if obj == "Sharpe":
                weights = np.matrix(w.value / k.value).T
            else:
                weights = np.matrix(w.value).T

            if self.sht == False:
                weights = np.abs(weights) / np.sum(np.abs(weights))

            for j in self.assetslist:
                portafolio[j].append(weights[0, self.assetslist.index(j)])

        except:
            pass

        optimum = pd.DataFrame(portafolio, index=["weights"], dtype=np.float64).T

        return optimum

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
            Risk measure used by the optimization model. The default is 'MV'.
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

        limits = pd.concat([w_min, w_max], axis=1)
        limits.columns = ["w_min", "w_max"]

        return limits

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
            Risk measure used by the optimization model. The default is 'MV'.
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
            mu = np.matrix(self.mu)
            sigma = np.matrix(self.cov)
            returns = np.matrix(self.returns)
            nav = np.matrix(self.nav)
        elif model == "FM":
            mu = np.matrix(self.mu_fm)
            if hist == False:
                sigma = np.matrix(self.cov_fm)
                returns = np.matrix(self.returns_fm)
                nav = np.matrix(self.nav_fm)
            elif hist == True:
                sigma = np.matrix(self.cov)
                returns = np.matrix(self.returns)
                nav = np.matrix(self.nav)
        elif model == "BL":
            mu = np.matrix(self.mu_bl)
            if hist == False:
                sigma = np.matrix(self.cov_bl)
            elif hist == True:
                sigma = np.matrix(self.cov)
            returns = np.matrix(self.returns)
            nav = np.matrix(self.nav)
        elif model == "BL_FM":
            mu = np.matrix(self.mu_bl_fm_2)
            if hist == False:
                sigma = np.matrix(self.cov_bl_fm_2)
                returns = np.matrix(self.returns_fm)
                nav = np.matrix(self.nav_fm)
            elif hist == True:
                sigma = np.matrix(self.cov)
                returns = np.matrix(self.returns)
                nav = np.matrix(self.nav)

        alpha1 = self.alpha

        limits = self.frontier_limits(model="Classic", rm=rm, rf=rf, hist=hist)

        w_min = np.matrix(limits.iloc[:, 0]).T
        w_max = np.matrix(limits.iloc[:, 1]).T

        ret_min = (mu * w_min).item()
        ret_max = (mu * w_max).item()

        if rm == "MV":
            risk_min = np.sqrt(w_min.T * sigma * w_min).item()
            risk_max = np.sqrt(w_max.T * sigma * w_max).item()
        elif rm == "MAD":
            risk_min = rk.MAD(returns * w_min)
            risk_max = rk.MAD(returns * w_max)
        elif rm == "MSV":
            risk_min = rk.SemiDeviation(returns * w_min)
            risk_max = rk.SemiDeviation(returns * w_max)
        elif rm == "CVaR":
            risk_min = rk.CVaR_Hist(returns * w_min, alpha1)
            risk_max = rk.CVaR_Hist(returns * w_max, alpha1)
        elif rm == "WR":
            risk_min = rk.WR(returns * w_min)
            risk_max = rk.WR(returns * w_max)
        elif rm == "FLPM":
            risk_min = rk.LPM(returns * w_min, rf, 1)
            risk_max = rk.LPM(returns * w_max, rf, 1)
        elif rm == "SLPM":
            risk_min = rk.LPM(returns * w_min, rf, 2)
            risk_max = rk.LPM(returns * w_max, rf, 2)
        elif rm == "MDD":
            risk_min = rk.MaxAbsDD(returns * w_min)
            risk_max = rk.MaxAbsDD(returns * w_max)
        elif rm == "ADD":
            risk_min = rk.AvgAbsDD(returns * w_min)
            risk_max = rk.AvgAbsDD(returns * w_max)
        elif rm == "CDaR":
            risk_min = rk.ConAbsDD(returns * w_min, alpha1)
            risk_max = rk.ConAbsDD(returns * w_max, alpha1)

        mus = np.linspace(ret_min, ret_max + (ret_max - ret_min) / (points), points + 1)

        risks = np.linspace(
            risk_min, risk_max + (risk_max - risk_min) / (points), points + 1
        )

        risk_lims = [
            "upperdev",
            "uppermad",
            "uppersdev",
            "upperCVaR",
            "upperwr",
            "upperflpm",
            "upperslpm",
            "uppermdd",
            "upperadd",
            "upperCDaR",
        ]

        risk_names = [
            "MV",
            "MAD",
            "MSV",
            "CVaR",
            "WR",
            "FLPM",
            "SLPM",
            "MDD",
            "ADD",
            "CDaR",
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
                n += 1
                frontier.append(w)
            except:
                pass

        setattr(self, risk_lims[item], None)
        frontier = pd.concat(frontier, axis=1)
        frontier.columns = list(range(len(risks)))

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
            "upperwr",
            "upperflpm",
            "upperslpm",
            "uppermdd",
            "upperadd",
            "upperCDaR",
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
