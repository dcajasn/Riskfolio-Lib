import sys
import numpy as np
import pandas as pd
import scipy.cluster.hierarchy as hr
from scipy.spatial.distance import squareform
import riskfolio.RiskFunctions as rk
import riskfolio.AuxFunctions as af
import riskfolio.Portfolio as pf


class HCPortfolio(object):
    r"""
    Class that creates a portfolio object with all properties needed to
    calculate optimal portfolios.

    Parameters
    ----------
    returns : DataFrame, optional
        A dataframe that containts the returns of the assets.
        The default is None.
    alpha : float, optional
        Significance level of CVaR, EVaR, CDaR and EDaR. The default is 0.05.
    """

    def __init__(self, returns=None, alpha=0.05):
        self._returns = returns
        self.alpha = alpha
        self.asset_order = None
        self.clusters = None
        self.cov = None
        self.corr = None
        self.corr_sorted = None

        # Solver params

        self.solvers = ["ECOS", "SCS", "OSQP", "CVXOPT"]
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

    # get naive-risk weights
    def _naive_risk(self, returns, cov, rm="MV", rf=0):
        assets = returns.columns.tolist()
        n = len(assets)

        if rm == "equal":
            weight = np.ones((n, 1)) * 1 / n
        else:
            inv_risk = np.zeros((n, 1))
            for i in assets:
                k = assets.index(i)
                w = np.zeros((n, 1))
                w[k, 0] = 1
                w = pd.DataFrame(w, columns=["weights"], index=assets)
                if rm == "vol":
                    risk = rk.Sharpe_Risk(
                        w, cov=cov, returns=returns, rm="MV", rf=rf, alpha=self.alpha
                    )
                else:
                    risk = rk.Sharpe_Risk(
                        w, cov=cov, returns=returns, rm=rm, rf=rf, alpha=self.alpha
                    )
                inv_risk[k, 0] = risk

            if rm == "MV":
                inv_risk = 1 / np.power(inv_risk, 2)
            else:
                inv_risk = 1 / inv_risk
            weight = inv_risk * (1 / np.sum(inv_risk))

        weight = weight.reshape(-1, 1)

        return weight

    # create hierarchical clustering
    def _hierarchical_clustering_hrp(self, linkage="single", leaf_order=True):

        # hierarchcial clustering
        dist = np.sqrt(
            np.clip((1.0 - self.corr) / 2.0, a_min=0.0, a_max=1.0)
        ).to_numpy()
        p_dist = squareform(dist, checks=False)
        clusters = hr.linkage(p_dist, method=linkage, optimal_ordering=leaf_order)

        return clusters

    # create hierarchical clustering
    def _hierarchical_clustering_herc(self, linkage="ward", max_k=10, leaf_order=True):

        # hierarchcial clustering
        dist = np.sqrt((1 - self.corr).round(8) / 2)
        dist = pd.DataFrame(dist, columns=self.corr.columns, index=self.corr.index)
        p_dist = squareform(dist, checks=False)
        clustering = hr.linkage(p_dist, method=linkage, optimal_ordering=leaf_order)

        # optimal number of clusters
        k = af.two_diff_gap_stat(self.corr, dist, clustering, max_k)

        return clustering, k

    # sort clustered items by distance
    def _seriation(self, clusters):
        return hr.leaves_list(clusters)

    # compute HRP weight allocation through recursive bisection
    def _recursive_bisection(self, sort_order, rm="MV", rf=0):
        weight = pd.Series(1, index=sort_order)  # set initial weights to 1
        items = [sort_order]

        while len(items) > 0:  # loop while weights is under 100%
            items = [
                i[j:k]
                for i in items
                for j, k in (
                    (0, len(i) // 2),
                    (len(i) // 2, len(i)),
                )  # get cluster indi
                if len(i) > 1
            ]

            # allocate weight to left and right cluster
            for i in range(0, len(items), 2):
                left_cluster = items[i]
                right_cluster = items[i + 1]

                # Left cluster
                left_cov = self.cov.iloc[left_cluster, left_cluster]
                left_returns = self.returns.iloc[:, left_cluster]
                left_weight = self._naive_risk(left_returns, left_cov, rm=rm, rf=rf)

                if rm == "vol":
                    left_risk = rk.Sharpe_Risk(
                        left_weight,
                        cov=left_cov,
                        returns=left_returns,
                        rm="MV",
                        rf=rf,
                        alpha=self.alpha,
                    )
                else:
                    left_risk = rk.Sharpe_Risk(
                        left_weight,
                        cov=left_cov,
                        returns=left_returns,
                        rm=rm,
                        rf=rf,
                        alpha=self.alpha,
                    )
                    if rm == "MV":
                        left_risk = np.power(left_risk, 2)

                # Right cluster
                right_cov = self.cov.iloc[right_cluster, right_cluster]
                right_returns = self.returns.iloc[:, right_cluster]
                right_weight = self._naive_risk(right_returns, right_cov, rm=rm, rf=rf)

                if rm == "vol":
                    right_risk = rk.Sharpe_Risk(
                        right_weight,
                        cov=right_cov,
                        returns=right_returns,
                        rm="MV",
                        rf=rf,
                        alpha=self.alpha,
                    )
                else:
                    right_risk = rk.Sharpe_Risk(
                        right_weight,
                        cov=right_cov,
                        returns=right_returns,
                        rm=rm,
                        rf=rf,
                        alpha=self.alpha,
                    )
                    if rm == "MV":
                        right_risk = np.power(right_risk, 2)

                # Allocate weight to clusters
                alpha = 1 - left_risk / (left_risk + right_risk)
                weight[left_cluster] *= alpha  # weight 1
                weight[right_cluster] *= 1 - alpha  # weight 2

        weight.index = self.asset_order

        return weight

    # compute HRP weight allocation through cluster-based bisection
    def _hierarchical_recursive_bisection(self, Z, rm="MV", rf=0, linkage="ward"):

        # Transform linkage to tree and reverse order
        root, nodes = hr.to_tree(Z, rd=True)
        nodes = nodes[::-1]
        weight = pd.Series(1, index=self.cov.index)  # Set initial weights to 1

        clusters_inds = hr.fcluster(Z, self.k, criterion="maxclust")
        clusters = {i: [] for i in range(min(clusters_inds), max(clusters_inds) + 1)}
        for i, v in enumerate(clusters_inds):
            clusters[v].append(i)

        # Loop through k clusters
        for i in nodes[: self.k - 1]:
            if i.is_leaf() == False:  # skip leaf-nodes
                left = i.get_left().pre_order()  # lambda i: i.id) # get left cluster
                right = i.get_right().pre_order()  # lambda i: i.id) # get right cluster
                left_set = set(left)
                right_set = set(right)
                left_risk = 0
                right_risk = 0

                # Allocate weight to clusters
                if rm == "equal":
                    w_1 = 0.5

                else:
                    for j in clusters.keys():
                        if set(clusters[j]).issubset(left_set):
                            # Left cluster
                            left_cov = self.cov.iloc[clusters[j], clusters[j]]
                            left_returns = self.returns.iloc[:, clusters[j]]
                            left_weight = self._naive_risk(
                                left_returns, left_cov, rm=rm, rf=rf
                            )

                            if rm == "vol":
                                left_risk_ = rk.Sharpe_Risk(
                                    left_weight,
                                    cov=left_cov,
                                    returns=left_returns,
                                    rm="MV",
                                    rf=rf,
                                    alpha=self.alpha,
                                )
                            else:
                                left_risk_ = rk.Sharpe_Risk(
                                    left_weight,
                                    cov=left_cov,
                                    returns=left_returns,
                                    rm=rm,
                                    rf=rf,
                                    alpha=self.alpha,
                                )
                                if rm == "MV":
                                    left_risk_ = np.power(left_risk_, 2)

                            left_risk += left_risk_

                        if set(clusters[j]).issubset(right_set):
                            # Right cluster
                            right_cov = self.cov.iloc[clusters[j], clusters[j]]
                            right_returns = self.returns.iloc[:, clusters[j]]
                            right_weight = self._naive_risk(
                                right_returns, right_cov, rm=rm, rf=rf
                            )

                            if rm == "vol":
                                right_risk_ = rk.Sharpe_Risk(
                                    right_weight,
                                    cov=right_cov,
                                    returns=right_returns,
                                    rm="MV",
                                    rf=rf,
                                    alpha=self.alpha,
                                )
                            else:
                                right_risk_ = rk.Sharpe_Risk(
                                    right_weight,
                                    cov=right_cov,
                                    returns=right_returns,
                                    rm=rm,
                                    rf=rf,
                                    alpha=self.alpha,
                                )
                                if rm == "MV":
                                    right_risk_ = np.power(right_risk_, 2)

                            right_risk += right_risk_

                    w_1 = 1 - left_risk / (left_risk + right_risk)

                weight[left] *= w_1  # weight 1
                weight[right] *= 1 - w_1  # weight 2

        # Get constituents of k clusters
        clustered_assets = pd.Series(
            hr.cut_tree(Z, n_clusters=self.k).flatten(), index=self.cov.index
        )

        # Multiply within-cluster weight with inter-cluster weight
        for i in range(self.k):
            cluster = clustered_assets.loc[clustered_assets == i]
            cluster_cov = self.cov.loc[cluster.index, cluster.index]
            cluster_returns = self.returns.loc[:, cluster.index]
            cluster_weights = pd.Series(
                self._naive_risk(cluster_returns, cluster_cov, rm=rm, rf=rf).flatten(),
                index=cluster_cov.index,
            )
            weight.loc[cluster_weights.index] *= cluster_weights

        return weight

    # Allocate weights
    def optimization(
        self,
        model="HRP",
        correlation="pearson",
        rm="MV",
        rf=0,
        linkage="single",
        k=None,
        max_k=10,
        leaf_order=True,
    ):
        r"""
        This method calculates the optimal portfolio according to the
        optimization model selected by the user.
        
        Parameters
        ----------
        model : str can be {'HRP' or 'HERC'}
            The hierarchical cluster portfolio model used for optimize the
            portfolio. The default is 'HRP'. Posible values are:

            - 'HRP': Hierarchical Risk Parity.
            - 'HERC': Hierarchical Equal Risk Contribution.

        correlation : str can be {'pearson', 'spearman' or 'distance'}.
            The correlation matrix used for create the clusters.
            The default is 'pearson'. Posible values are:

            - 'pearson': pearson correlation matrix.
            - 'spearman': spearman correlation matrix.
            - 'abs_pearson': absolute value pearson correlation matrix.
            - 'abs_spearman': absolute value spearman correlation matrix.
            - 'distance': distance correlation matrix.

        rm : str, optional
            The risk measure used to optimze the portfolio.
            The default is 'MV'. Posible values are:
            
            - 'vol': Standard Deviation.
            - 'MV': Variance.
            - 'MAD': Mean Absolute Deviation.
            - 'MSV': Semi Standard Deviation.
            - 'FLPM': First Lower Partial Moment (Omega Ratio).
            - 'SLPM': Second Lower Partial Moment (Sortino Ratio).
            - 'VaR': Value at Risk.
            - 'CVaR': Conditional Value at Risk.
            - 'EVaR': Entropic Value at Risk.
            - 'WR': Worst Realization (Minimax)
            - 'MDD': Maximum Drawdown of uncompounded cumulative returns (Calmar Ratio).
            - 'ADD': Average Drawdown of uncompounded cumulative returns.
            - 'DaR': Drawdown at Risk of uncompounded cumulative returns.
            - 'CDaR': Conditional Drawdown at Risk of uncompounded cumulative returns.
            - 'EDaR': Entropic Drawdown at Risk of uncompounded cumulative returns.
            - 'UCI': Ulcer Index of uncompounded cumulative returns.
            - 'MDD_Rel': Maximum Drawdown of compounded cumulative returns (Calmar Ratio).
            - 'ADD_Rel': Average Drawdown of compounded cumulative returns.
            - 'DaR_Rel': Drawdown at Risk of compounded cumulative returns.
            - 'CDaR_Rel': Conditional Drawdown at Risk of compounded cumulative returns.
            - 'EDaR_Rel': Entropic Drawdown at Risk of compounded cumulative returns.
            - 'UCI_Rel': Ulcer Index of compounded cumulative returns.
                
        rf : float, optional
            Risk free rate, must be in the same period of assets returns.
            The default is 0.
        linkage : string, optional
            Linkage method of hierarchical clustering, see `linkage <https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html?highlight=linkage#scipy.cluster.hierarchy.linkage>`_ for more details.
            The default is 'single'. Posible values are:

            - 'single'.
            - 'complete'.
            - 'average'.
            - 'weighted'.
            - 'centroid'.
            - 'median'.
            - 'ward'.
        
        k : int, optional
            Number of clusters. This value is took instead of the optimal number
            of clusters calculated with the two difference gap statistic.
            The default is None.
        max_k : int, optional
            Max number of clusters used by the two difference gap statistic
            to find the optimal number of clusters. The default is 10.
        
        Returns
        -------
        w : DataFrame
            The weights of optimal portfolio.

        """

        # Correlation matrix from covariance matrix
        self.cov = self.returns.cov()
        if correlation in {"pearson", "spearman"}:
            self.corr = self.returns.corr(method=correlation)
        if correlation in {"abs_pearson", "abs_spearman"}:
            self.corr = np.abs(self.returns.corr(method=correlation[4:]))
        elif correlation == "distance":
            self.corr = af.dcorr_matrix(self.returns)

        # Step-1: Tree clustering
        if model == "HRP":
            self.clusters = self._hierarchical_clustering_hrp(
                linkage, leaf_order=leaf_order
            )
        elif model == "HERC":
            self.clusters, self.k = self._hierarchical_clustering_herc(
                linkage, max_k, leaf_order=leaf_order
            )
            if k is not None:
                self.k = int(k)

        # Step-2: Seriation (Quasi-Diagnalization)
        self.sort_order = self._seriation(self.clusters)
        asset_order = self.assetslist
        asset_order[:] = [self.assetslist[i] for i in self.sort_order]
        self.asset_order = asset_order
        self.corr_sorted = self.corr.reindex(
            index=self.asset_order, columns=self.asset_order
        )

        # Step-3: Recursive bisection
        if model == "HRP":
            weights = self._recursive_bisection(self.sort_order, rm=rm, rf=rf)
        elif model == "HERC":
            weights = self._hierarchical_recursive_bisection(
                self.clusters, rm=rm, rf=rf, linkage=linkage
            )

        weights = weights.loc[self.assetslist].to_frame()

        return weights
