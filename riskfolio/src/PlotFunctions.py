""""""  #

"""
Copyright (c) 2020-2025, Dany Cajas
All rights reserved.
This work is licensed under BSD 3-Clause "New" or "Revised" License.
License available at https://github.com/dcajasn/Riskfolio-Lib/blob/master/LICENSE.txt
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.lines as mlines
import matplotlib.ticker as mticker
from matplotlib import cm, colors
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
import scipy.stats as st
import scipy.cluster.hierarchy as hr
from scipy.spatial.distance import squareform
import networkx as nx
import riskfolio.src.RiskFunctions as rk
import riskfolio.src.AuxFunctions as af
import riskfolio.src.DBHT as db
import riskfolio.src.ConstraintsFunctions as ct


__all__ = [
    "plot_series",
    "plot_frontier",
    "plot_pie",
    "plot_bar",
    "plot_frontier_area",
    "plot_risk_con",
    "plot_factor_risk_con",
    "plot_hist",
    "plot_range",
    "plot_drawdown",
    "plot_table",
    "plot_clusters",
    "plot_dendrogram",
    "plot_network",
    "plot_network_allocation",
    "plot_clusters_network",
    "plot_clusters_network_allocation",
    "plot_BrinsonAttribution",
]

rm_names = [
    "Standard Deviation",
    "Square Root Kurtosis",
    "Mean Absolute Deviation",
    "Gini Mean Difference",
    "Semi Standard Deviation",
    "Square Root Semi Kurtosis",
    "First Lower Partial Moment",
    "Second Lower Partial Moment",
    "Value at Risk",
    "Conditional Value at Risk",
    "Tail Gini",
    "Entropic Value at Risk",
    "Relativistic Value at Risk",
    "Worst Realization",
    "Conditional Value at Risk Range",
    "Tail Gini Range",
    "Entropic Value at Risk Range",
    "Relativistic Value at Risk Range",
    "Range",
    "Max Drawdown",
    "Average Drawdown",
    "Drawdown at Risk",
    "Conditional Drawdown at Risk",
    "Entropic Drawdown at Risk",
    "Relativistic Drawdown at Risk",
    "Ulcer Index",
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
    "VaR",
    "CVaR",
    "TG",
    "EVaR",
    "RLVaR",
    "WR",
    "CVRG",
    "TGRG",
    "EVRG",
    "RVRG",
    "RG",
    "MDD",
    "ADD",
    "DaR",
    "CDaR",
    "EDaR",
    "RLDaR",
    "UCI",
]


def plot_series(returns, w, cmap="tab20", n_colors=20, height=6, width=10, ax=None):
    r"""
    Create a chart with the compounded cumulative of the portfolios.

    Parameters
    ----------
    returns : DataFrame of shape (n_samples, n_assets)
        Assets returns DataFrame, where n_samples is the number of
        observations and n_assets is the number of assets.
    w : DataFrame or Series of shape (n_assets, 1)
        Portfolio weights, where n_assets is the number of assets.
    cmap : cmap, optional
        Colorscale used to plot each portfolio compounded cumulative return.
        The default is 'tab20'.
    n_colors : int, optional
        Number of distinct colors per color cycle. If there are more assets
        than n_colors, the chart is going to start to repeat the color cycle.
        The default is 20.
    height : float, optional
        Height of the image in inches. The default is 6.
    width : float, optional
        Width of the image in inches. The default is 10.
    ax : matplotlib axis, optional
        If provided, plot on this axis. The default is None.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    Returns
    -------
    ax : matplotlib axis
        Returns the Axes object with the plot for further tweaking.

    Example
    -------
    ::

        ax = rp.plot_series(returns=Y,
                            w=ws,
                            cmap='tab20',
                            height=6,
                            width=10,
                            ax=None)

    .. image:: images/Port_Series.png

    """

    if not isinstance(returns, pd.DataFrame):
        raise ValueError("returns must be a DataFrame")

    if not isinstance(w, pd.DataFrame):
        if isinstance(w, pd.Series):
            w_ = w.to_frame()
        else:
            raise ValueError("w must be a DataFrame or Series.")
    else:
        w_ = w.copy()

    if returns.columns.tolist() != w_.index.tolist():
        if returns.columns.tolist() == w_.index.tolist():
            w_ = w_.T
        else:
            raise ValueError("returns and w must have the same assets.")

    if ax is None:
        fig = plt.gcf()
        ax = fig.gca()
        fig.set_figwidth(width)
        fig.set_figheight(height)
    else:
        fig = ax.get_figure()

    ax.grid(linestyle=":")
    title = "Historical Compounded Cumulative Returns"
    ax.set_title(title)

    labels = w_.columns.tolist()
    index = returns.index.tolist()

    colormap = cm.get_cmap(cmap)
    colormap = colormap(np.linspace(0, 1, int(n_colors)))

    if cmap == "gist_rainbow":
        colormap = colormap[::-1]

    cycle = plt.cycler("color", colormap)
    ax.set_prop_cycle(cycle)

    for i in range(len(labels)):
        a = np.array(returns, ndmin=2) @ np.array(w_[labels[i]], ndmin=2).T
        prices = 1 + np.insert(a, 0, 0, axis=0)
        prices = np.cumprod(prices, axis=0)
        prices = np.ravel(prices).tolist()
        del prices[0]

        ax.plot_date(index, prices, "-", label=labels[i])

    ax.xaxis.set_major_locator(mdates.AutoDateLocator(tz=None, minticks=5, maxticks=10))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))

    ticks_loc = ax.get_yticks().tolist()
    ax.set_yticks(ax.get_yticks().tolist())
    ax.set_yticklabels(["{:3.2f}".format(x) for x in ticks_loc])
    ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))

    try:
        fig.set_layout_engine(layout="constrained")
    except:
        pass

    return ax


def plot_frontier(
    w_frontier,
    returns,
    mu=None,
    cov=None,
    rm="MV",
    kelly=False,
    rf=0,
    alpha=0.05,
    a_sim=100,
    beta=None,
    b_sim=None,
    kappa=0.30,
    kappa_g=None,
    solver="CLARABEL",
    cmap="viridis",
    w=None,
    label="Portfolio",
    marker="*",
    s=16,
    c="r",
    height=6,
    width=10,
    t_factor=252,
    ax=None,
):
    r"""
    Creates a plot of the efficient frontier for a risk measure specified by
    the user.

    Parameters
    ----------
    w_frontier : DataFrame
        Portfolio weights of some points in the efficient frontier.
    mu : DataFrame of shape (1, n_assets)
        Vector of expected returns, where n_assets is the number of assets.
    cov : DataFrame of shape (n_assets, n_assets)
        Covariance matrix, where n_assets is the number of assets.
    returns : DataFrame of shape (n_samples, n_assets)
        Assets returns DataFrame, where n_samples is the number of
        observations and n_assets is the number of assets.
    rm : str, optional
        The risk measure used to estimate the frontier.
        The default is 'MV'. Possible values are:

        - 'MV': Standard Deviation.
        - 'KT': Square Root Kurtosis.
        - 'MAD': Mean Absolute Deviation.
        - 'MSV': Semi Standard Deviation.
        - 'SKT': Square Root Semi Kurtosis.
        - 'FLPM': First Lower Partial Moment (Omega Ratio).
        - 'SLPM': Second Lower Partial Moment (Sortino Ratio).
        - 'CVaR': Conditional Value at Risk.
        - 'TG': Tail Gini.
        - 'EVaR': Entropic Value at Risk.
        - 'RLVaR': Relativistic Value at Risk.
        - 'WR': Worst Realization (Minimax).
        - 'CVRG': CVaR range of returns.
        - 'TGRG': Tail Gini range of returns.
        - 'EVRG': EVaR range of returns.
        - 'RVRG': RLVaR range of returns.
        - 'RG': Range of returns.
        - 'MDD': Maximum Drawdown of uncompounded returns (Calmar Ratio).
        - 'ADD': Average Drawdown of uncompounded cumulative returns.
        - 'DaR': Drawdown at Risk of uncompounded cumulative returns.
        - 'CDaR': Conditional Drawdown at Risk of uncompounded cumulative returns.
        - 'EDaR': Entropic Drawdown at Risk of uncompounded cumulative returns.
        - 'RLDaR': Relativistic Drawdown at Risk of uncompounded cumulative returns.
        - 'UCI': Ulcer Index of uncompounded cumulative returns.

    kelly : bool, optional
        Method used to calculate mean return. Possible values are False for
        arithmetic mean return and True for mean logarithmic return. The default
        is False.
    rf : float, optional
        Risk free rate or minimum acceptable return. The default is 0.
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
        Deformation parameter of RLVaR and RLDaR for losses, must be between 0 and 1. The default is 0.30.
    kappa_g : float, optional
        Deformation parameter of RLVaR for gains, must be between 0 and 1.
        The default is None.
    solver: str, optional
        Solver available for CVXPY that supports power cone programming. Used to calculate RLVaR and RLDaR.
        The default value is 'CLARABEL'.
    cmap : cmap, optional
        Colorscale that represents the risk adjusted return ratio.
        The default is 'viridis'.
    w : DataFrame  or Series of shape (n_assets, 1), optional
        A portfolio specified by the user. The default is None.
    label : str or list, optional
        Name or list of names of portfolios that appear on plot legend.
        The default is 'Portfolio'.
    marker : str, optional
        Marker of w. The default is "*".
    s : float, optional
        Size of marker. The default is 16.
    c : str, optional
        Color of marker. The default is 'r'.
    height : float, optional
        Height of the image in inches. The default is 6.
    width : float, optional
        Width of the image in inches. The default is 10.
    t_factor : float, optional
        Factor used to annualize expected return and expected risks for
        risk measures based on returns (not drawdowns). The default is 252.

        .. math::

            \begin{align}
            \text{Annualized Return} & = \text{Return} \, \times \, \text{t_factor} \\
            \text{Annualized Risk} & = \text{Risk} \, \times \, \sqrt{\text{t_factor}}
            \end{align}

    ax : matplotlib axis, optional
        If provided, plot on this axis. The default is None.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    Returns
    -------
    ax : matplotlib Axes
        Returns the Axes object with the plot for further tweaking.

    Example
    -------
    ::

        label = 'Max Risk Adjusted Return Portfolio'
        mu = port.mu
        cov = port.cov
        returns = port.returns

        ax = rp.plot_frontier(w_frontier=ws,
                              mu=mu,
                              cov=cov,
                              returns=Y,
                              rm=rm,
                              rf=0,
                              alpha=0.05,
                              cmap='viridis',
                              w=w1,
                              label=label,
                              marker='*',
                              s=16,
                              c='r',
                              height=6,
                              width=10,
                              t_factor=252,
                              ax=None)

    .. image:: images/MSV_Frontier.png


    """

    if not isinstance(w_frontier, pd.DataFrame):
        raise ValueError("w_frontier must be a DataFrame.")

    if not isinstance(returns, pd.DataFrame):
        raise ValueError("returns must be a DataFrame.")
    else:
        if returns.columns.tolist() != w_frontier.index.tolist():
            if returns.columns.tolist() != w_frontier.columns.tolist():
                raise ValueError("returns and w_frontier must have same assets.")
            else:
                w_frontier_ = w_frontier.T.copy()
        else:
            w_frontier_ = w_frontier.copy()

    if mu is None:
        mu_ = returns.mean()
    elif isinstance(mu, pd.DataFrame):
        if mu.shape[0] == 1 or mu.shape[1] == 1:
            mu_ = mu.squeeze()
        else:
            raise ValueError("mu must be a DataFrame or Series of one dimension.")
    elif isinstance(mu, pd.Series):
        mu_ = mu.copy()
    else:
        raise ValueError("mu must be a DataFrame or Series.")

    if cov is None:
        cov = returns.cov()
    elif isinstance(cov, pd.DataFrame):
        if cov.index.tolist() != cov.columns.tolist():
            raise ValueError(
                "cov must be a square DataFrame with samen labels in indexes and columns."
            )
        elif cov.index.tolist() != w_frontier_.index.tolist():
            raise ValueError("cov and w_frontier must have the same assets.")
    else:
        raise ValueError("cov must be a square DataFrame.")

    if w is not None:
        if not isinstance(w, pd.DataFrame):
            if isinstance(w, pd.Series):
                w_ = w.to_frame()
            else:
                raise ValueError("w must be a DataFrame or Series")
        else:
            if returns.columns.tolist() == w.columns.tolist():
                w_ = w.T.copy()
            elif returns.columns.tolist() == w.index.tolist():
                w_ = w.copy()
            else:
                raise ValueError("returns and w must have the same assets.")

    if beta is None:
        beta = alpha
    if b_sim is None:
        b_sim = a_sim
    if kappa_g is None:
        kappa_g = kappa

    width_ratios = [0.97, 0.03]

    if ax is None:
        fig = plt.gcf()
        ax = fig.gca()
        ax.axis("off")
        gs = GridSpec(nrows=1, ncols=2, figure=fig, width_ratios=width_ratios)
        axes = []
        axes.append(fig.add_subplot(gs[0]))
        axes.append(fig.add_subplot(gs[1]))
        fig.set_figwidth(width)
        fig.set_figheight(height)
    else:
        if isinstance(ax, plt.Axes):
            ax.axis("off")
            fig = ax.get_figure()
            if hasattr(ax, "get_subplotspec"):
                subplot_spec = ax.get_subplotspec()
                gs0 = subplot_spec.get_gridspec()
                gs = GridSpecFromSubplotSpec(
                    nrows=1, ncols=2, width_ratios=width_ratios, subplot_spec=gs0[0]
                )
            else:
                gs = GridSpec(nrows=1, ncols=2, figure=fig, width_ratios=width_ratios)
            axes = []
            axes.append(fig.add_subplot(gs[0]))
            axes.append(fig.add_subplot(gs[1]))
        else:
            raise TypeError("ax must be a matplotlib axes object.")

    mu_ = np.array(mu_, ndmin=2)

    ax0 = axes[0]
    if kelly == False:
        ax0.set_ylabel("Expected Arithmetic Return")
    elif kelly == True:
        ax0.set_ylabel("Expected Logarithmic Return")

    item = rmeasures.index(rm)
    if rm in [
        "CVaR",
        "TG",
        "EVaR",
        "RLVaR",
        "CVRG",
        "TGRG",
        "EVRG",
        "RVRG",
        "CDaR",
        "EDaR",
        "RLDaR",
    ]:
        x_label = (
            rm_names[item] + " (" + rm + ")" + " $\\alpha = $" + "{0:.2%}".format(alpha)
        )
    else:
        x_label = rm_names[item] + " (" + rm + ")"
    if rm in ["CVRG", "TGRG", "EVRG", "RVRG"]:
        x_label += ", $\\beta = $" + "{0:.2%}".format(beta)
    if rm in ["RLVaR", "RLDaR", "RVRG"]:
        x_label += ", $\\kappa = $" + "{0:.2}".format(kappa)
    if rm in ["RLVaR", "RLDaR", "RVRG"]:
        x_label += ", $\\kappa_g = $" + "{0:.2}".format(kappa_g)

    ax0.set_xlabel("Expected Risk - " + x_label)

    title = "Efficient Frontier Mean - " + x_label
    ax0.set_title(title)

    X1 = []
    Y1 = []
    Z1 = []

    for i in range(w_frontier_.shape[1]):
        try:
            weights = np.array(w_frontier_.iloc[:, i], ndmin=2).T
            risk = rk.Sharpe_Risk(
                returns=returns,
                w=weights,
                cov=cov,
                rm=rm,
                rf=rf,
                alpha=alpha,
                a_sim=a_sim,
                beta=beta,
                b_sim=b_sim,
                kappa=kappa,
                solver=solver,
            )

            if kelly == False:
                ret = mu_ @ weights
            elif kelly == True:
                ret = 1 / returns.shape[0] * np.sum(np.log(1 + returns @ weights))
            ret = ret.item() * t_factor

            if rm not in ["MDD", "ADD", "CDaR", "EDaR", "RLDaR", "UCI"]:
                risk = risk * t_factor**0.5

            ratio = (ret - rf) / risk

            X1.append(risk)
            Y1.append(ret)
            Z1.append(ratio)
        except:
            pass

    ax_scatter = ax0.scatter(X1, Y1, c=Z1, cmap=cmap)

    if w is not None:
        if isinstance(label, str):
            label = [label]

        if label is None:
            label = w_.columns.tolist()

        if w_.shape[1] != len(label):
            label = w_.columns.tolist()

        label = [
            v + " " + str(label[:i].count(v) + 1) if label.count(v) > 1 else v
            for i, v in enumerate(label)
        ]

        if isinstance(c, str):
            colormap = np.array(colors.to_rgba(c)).reshape(1, -1)
        elif c is None:
            colormap = np.array(colors.to_rgba("red")).reshape(1, -1)

        elif isinstance(c, list):
            colormap = [list(colors.to_rgba(i)) for i in c]
            colormap = np.array(colormap)

        if len(label) != colormap.shape[0]:
            colormap = cm.get_cmap("tab20")
            colormap = colormap(np.linspace(0, 1, 20))
            colormap = np.vstack(
                [colormap[6:8], colormap[2:6], colormap[8:], colormap[0:2]]
            )

        n_repeats = int(len(label) // 20 + 1)
        if n_repeats > 1:
            colormap = np.vstack([colormap] * n_repeats)

        for i in range(w_.shape[1]):
            weights = w_.iloc[:, i].to_numpy().reshape(-1, 1)
            risk = rk.Sharpe_Risk(
                returns=returns,
                w=weights,
                cov=cov,
                rm=rm,
                rf=rf,
                alpha=alpha,
                a_sim=a_sim,
                beta=beta,
                b_sim=b_sim,
                kappa=kappa,
                solver=solver,
            )
            if kelly == False:
                ret = mu_ @ weights
            elif kelly == True:
                ret = 1 / returns.shape[0] * np.sum(np.log(1 + returns @ weights))
            ret = ret.item() * t_factor

            if rm not in ["MDD", "ADD", "CDaR", "EDaR", "RLDaR", "UCI"]:
                risk = risk * t_factor**0.5

            color = colormap[i].reshape(1, -1)
            ax0.scatter(risk, ret, marker=marker, s=s**2, c=color, label=label[i])

        ax0.legend(loc="upper left")

    xmin = np.min(X1) - np.abs(np.max(X1) - np.min(X1)) * 0.1
    xmax = np.max(X1) + np.abs(np.max(X1) - np.min(X1)) * 0.1
    ymin = np.min(Y1) - np.abs(np.max(Y1) - np.min(Y1)) * 0.1
    ymax = np.max(Y1) + np.abs(np.max(Y1) - np.min(Y1)) * 0.1

    ax0.set_ylim(ymin, ymax)
    ax0.set_xlim(xmin, xmax)

    ax0.xaxis.set_major_locator(plt.AutoLocator())

    ticks_loc = ax0.get_yticks().tolist()
    ax0.set_yticks(ax0.get_yticks().tolist())
    ax0.set_yticklabels(["{:.2%}".format(x) for x in ticks_loc])
    ticks_loc = ax0.get_xticks().tolist()
    ax0.set_xticks(ax0.get_xticks().tolist())
    ax0.set_xticklabels(["{:.2%}".format(x) for x in ticks_loc])

    ax0.tick_params(axis="y", direction="in")
    ax0.tick_params(axis="x", direction="in")

    ax0.grid(linestyle=":")

    axcolor = axes[1]
    colorbar = fig.colorbar(ax_scatter, cax=axcolor)
    colorbar.set_label("Risk Adjusted Return Ratio")

    try:
        fig.set_layout_engine(layout="constrained")
    except:
        pass

    return ax


def plot_pie(
    w,
    title="",
    others=0.05,
    nrow=25,
    cmap="tab20",
    n_colors=20,
    height=6,
    width=8,
    ax=None,
):
    r"""
    Create a pie chart with portfolio weights.

    Parameters
    ----------
    w : DataFrame or Series of shape (n_assets, 1)
        Portfolio weights, where n_assets is the number of assets.
    title : str, optional
        Title of the chart. The default is "".
    others : float, optional
        Percentage of others section. The default is 0.05.
    nrow : int, optional
        Number of rows of the legend. The default is 25.
    cmap : cmap, optional
        Color scale used to plot each asset weight.
        The default is 'tab20'.
    n_colors : int, optional
        Number of distinct colors per color cycle. If there are more assets
        than n_colors, the chart is going to start to repeat the color cycle.
        The default is 20.
    height : float, optional
        Height of the image in inches. The default is 10.
    width : float, optional
        Width of the image in inches. The default is 10.
    ax : matplotlib axis, optional
        If provided, plot on this axis. The default is None.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    Returns
    -------
    ax :  matplotlib axis.
        Returns the Axes object with the plot for further tweaking.

    Example
    -------
    ::

        ax = rp.plot_pie(w=w1,
                         title='Portfolio',
                         height=6,
                         width=10,
                         cmap="tab20",
                         ax=None)

    .. image:: images/Pie_Chart.png


    """

    if not isinstance(w, pd.DataFrame):
        if isinstance(w, pd.Series):
            w_ = w.to_frame()
        else:
            raise ValueError("w must be a one column DataFrame or Series")
    else:
        if w.shape[0] == 1:
            w_ = w.T.copy()
        elif w.shape[1] == 1:
            w_ = w.copy()
        else:
            raise ValueError("w must be a one column DataFrame or Series")

    if ax is None:
        ax = plt.gca()
        fig = plt.gcf()
        fig.set_figwidth(width)
        fig.set_figheight(height)
    else:
        fig = ax.get_figure()

    labels = w_.index.tolist()
    sizes = w_.iloc[:, 0].tolist()
    abs_sizes = [np.absolute(s) for s in sizes]
    sizes2 = pd.DataFrame([labels, abs_sizes, sizes]).T
    sizes2.columns = ["labels", "abs_values", "values"]
    sizes2 = sizes2.sort_values(by=["abs_values"], ascending=False)
    sizes2.index = [i for i in range(0, len(labels))]
    sizes3 = sizes2.cumsum()
    sizes3["abs_values"] = sizes3["abs_values"] / sizes3["abs_values"].max()
    l = sizes3[sizes3["abs_values"] >= 1 - others].index.tolist()[0]

    if l >= 0:
        a1 = sizes2["abs_values"].sum() - sizes2[sizes2.index <= l]["abs_values"].sum()
        a2 = sizes2["values"].sum() - sizes2[sizes2.index <= l]["values"].sum()
        item = pd.DataFrame(["Others", a1, a2]).T
        item.columns = ["labels", "abs_values", "values"]
        sizes2 = sizes2[sizes2.index <= l]
        sizes2 = pd.concat([sizes2, item], axis=0)

    abs_sizes = sizes2["abs_values"].tolist()
    sizes = sizes2["values"].tolist()
    labels = sizes2["labels"].tolist()
    sizes2 = ["{0:.1%}".format(i) for i in sizes]

    if title == "":
        title = "Portfolio Composition"

    limit = np.round(np.min(sizes), 4)
    if limit < 0:
        title += " (Areas in Absolute Values)"

    ax.set_title(title)

    colormap = cm.get_cmap(cmap)
    colormap = colormap(np.linspace(0, 1, n_colors))

    if cmap == "gist_rainbow":
        colormap = colormap[::-1]

    cycle = plt.cycler("color", colormap)
    ax.set_prop_cycle(cycle)

    size = 0.4

    # set up style cycles

    wedges, texts = ax.pie(
        abs_sizes,
        radius=1,
        wedgeprops=dict(width=size, edgecolor="black"),
        startangle=-15,
        normalize=True,
    )

    # Equal aspect ratio ensures that pie is drawn as a circle.

    ax.axis("equal")

    n = int(np.ceil(l / nrow))
    if n == 0:
        n += 1

    ax.legend(wedges, labels, loc="center left", bbox_to_anchor=(1, 0.5), ncol=n)

    bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72)
    kw = dict(
        xycoords="data",
        textcoords="data",
        arrowprops=dict(arrowstyle="-"),
        bbox=bbox_props,
        zorder=0,
        va="center",
    )

    for i, p in enumerate(wedges):
        ang = (p.theta2 - p.theta1) / 2.0 + p.theta1
        y = np.sin(np.deg2rad(ang))
        x = np.cos(np.deg2rad(ang))
        horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
        connectionstyle = "angle,angleA=0,angleB={}".format(ang)
        kw["arrowprops"].update({"connectionstyle": connectionstyle})
        name = str(labels[i]) + " " + str(sizes2[i])
        ax.annotate(
            name,
            xy=(x, y),
            xytext=(1.1 * np.sign(x), 1.1 * y),
            horizontalalignment=horizontalalignment,
            **kw
        )

    try:
        fig.set_layout_engine(layout="constrained")
    except:
        pass

    return ax


def plot_bar(
    w,
    title="",
    kind="v",
    others=0.05,
    nrow=25,
    cpos="tab:green",
    cneg="darkorange",
    cothers="dodgerblue",
    height=6,
    width=10,
    ax=None,
):
    r"""
    Create a bar chart with portfolio weights.

    Parameters
    ----------
    w : DataFrame or Series of shape (n_assets, 1)
        Portfolio weights, where n_assets is the number of assets.
    title : str, optional
        Title of the chart. The default is "".
    kind : str, optional
        Kind of bar plot, "v" for vertical bars and "h" for horizontal bars.
        The default is "v".
    others : float, optional
        Percentage of others section. The default is 0.05.
    nrow : int, optional
        Max number of bars that be plotted. The default is 25.
    cpos : str, optional
        Color for positives weights. The default is 'tab:green'.
    cneg : str, optional
        Color for negatives weights. The default is 'darkorange'.
    cothers : str, optional
        Color for others bar. The default is 'dodgerblue'.
    height : float, optional
        Height of the image in inches. The default is 10.
    width : float, optional
        Width of the image in inches. The default is 10.
    ax : matplotlib axis, optional
        If provided, plot on this axis. The default is None.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    Returns
    -------
    ax :  matplotlib axis.
        Returns the Axes object with the plot for further tweaking.

    Example
    -------
    ::

        ax = rp.plot_bar(w1,
                         title='Portfolio',
                         kind="v",
                         others=0.05,
                         nrow=25,
                         height=6,
                         width=10,
                         ax=None)

    .. image:: images/Bar_Chart.png


    """

    if not isinstance(w, pd.DataFrame):
        if isinstance(w, pd.Series):
            w_ = w.to_frame()
        else:
            raise ValueError("w must be a one column DataFrame or Series")
    else:
        if w.shape[0] == 1:
            w_ = w.T.copy()
        elif w.shape[1] == 1:
            w_ = w.copy()
        else:
            raise ValueError("w must be a one column DataFrame or Series")

    if ax is None:
        fig = plt.gcf()
        ax = fig.gca()
        fig.set_figwidth(width)
        fig.set_figheight(height)
    else:
        fig = ax.get_figure()

    labels = w_.index.tolist()
    sizes = w_.to_numpy().flatten().tolist()
    abs_sizes = [np.absolute(s) for s in sizes]
    sizes2 = pd.DataFrame([labels, abs_sizes, sizes]).T
    sizes2.columns = ["labels", "abs_values", "values"]
    sizes2 = sizes2.sort_values(by=["abs_values"], ascending=False)
    sizes2.index = [i for i in range(0, len(labels))]
    sizes3 = sizes2.cumsum()
    sizes3["abs_values"] = sizes3["abs_values"] / sizes3["abs_values"].max()

    l1 = sizes3[sizes3["abs_values"] >= 1 - others].index.tolist()
    if len(l1) > 0:
        l1 = l1[0]
    else:
        l1 = -1

    l2 = sizes2[sizes2["abs_values"] < 0.01].index.tolist()
    if len(l2) > 0:
        l2 = l2[0]
    else:
        l2 = -1

    if l1 > int(nrow):
        a1 = sizes2["abs_values"].sum() - sizes2[sizes2.index <= l1]["abs_values"].sum()
        a2 = sizes2["values"].sum() - sizes2[sizes2.index <= l1]["values"].sum()
        item = pd.DataFrame(["Others", a1, a2]).T
        item.columns = ["labels", "abs_values", "values"]
        sizes2 = sizes2[sizes2.index <= l1]
        sizes2 = sizes2.sort_values(by=["values"], ascending=False)
        sizes2 = pd.concat([sizes2, item], axis=0)
    elif l2 > 0:
        a1 = sizes2["abs_values"].sum() - sizes2[sizes2.index <= l2]["abs_values"].sum()
        a2 = sizes2["values"].sum() - sizes2[sizes2.index <= l2]["values"].sum()
        item = pd.DataFrame(["Others", a1, a2]).T
        item.columns = ["labels", "abs_values", "values"]
        sizes2 = sizes2[sizes2.index <= l2]
        sizes2 = sizes2.sort_values(by=["values"], ascending=False)
        sizes2 = pd.concat([sizes2, item], axis=0)
    else:
        sizes2 = sizes2.sort_values(by=["values"], ascending=False)

    sizes = sizes2["values"].tolist()
    labels = sizes2["labels"].tolist()
    sizes2 = ["{0:.1%}".format(i) for i in sizes]

    if title == "":
        title = "Portfolio Composition"

    ax.set_title(title)

    if kind == "v":
        sizes = np.array(sizes)
        labels = np.array(labels)
        ax.bar(labels, np.where(sizes >= 0, sizes, 0), color=cpos, width=0.5)
        ax.bar(labels, np.where(sizes < 0, sizes, 0), color=cneg, width=0.5)

        if l1 > int(nrow):
            ax.bar(
                labels, np.where(labels == "Others", sizes, 0), color=cothers, width=0.5
            )
            b = "Others (Sum Abs < " + "{:.1%}".format(others) + ")"
        elif l2 > 0:
            ax.bar(
                labels, np.where(labels == "Others", sizes, 0), color=cothers, width=0.5
            )
            b = "Others (Abs < " + "{:.1%}".format(0.01) + ")"
        else:
            b = None

        ticks_loc = ax.get_yticks().tolist()
        ax.set_yticks(ax.get_yticks().tolist())
        ax.set_yticklabels(["{:.2%}".format(x) for x in ticks_loc])
        ax.set_xlim(-0.5, len(sizes) - 0.5)

        r = plt.gcf().canvas.get_renderer()
        transf = ax.transData.inverted()

        for i, v in enumerate(sizes):
            t = ax.text(i, v, sizes2[i], color="black")
            bb = t.get_window_extent(renderer=r)
            bb = bb.transformed(transf)
            h_text = bb.height
            x_text = bb.x0 - (bb.x1 - bb.x0) * 0.4
            y_text = bb.y0
            if v >= 0:
                t.set_position((x_text, y_text + h_text * 0.8))
            else:
                t.set_position((x_text, y_text - h_text))

    elif kind == "h":
        sizes.reverse()
        labels.reverse()
        sizes2.reverse()
        sizes = np.array(sizes)
        labels = np.array(labels)
        ax.barh(labels, np.where(sizes >= 0, sizes, 0), color=cpos, height=0.5)
        ax.barh(labels, np.where(sizes < 0, sizes, 0), color=cneg, height=0.5)

        if l1 > int(nrow):
            ax.barh(
                labels,
                np.where(labels == "Others", sizes, 0),
                color=cothers,
                height=0.5,
            )
            b = "Others (Sum Abs < " + "{:.1%}".format(others) + ")"
        elif l2 > 0:
            ax.barh(
                labels,
                np.where(labels == "Others", sizes, 0),
                color=cothers,
                height=0.5,
            )
            b = "Others (Abs < " + "{:.1%}".format(0.01) + ")"
        else:
            b = None

        ticks_loc = ax.get_xticks().tolist()
        ax.set_xticks(ax.get_xticks().tolist())
        ax.set_xticklabels(["{:.2%}".format(x) for x in ticks_loc])
        ax.set_ylim(-0.5, len(sizes) - 0.5)

        r = plt.gcf().canvas.get_renderer()
        transf = ax.transData.inverted()

        for i, v in enumerate(sizes):
            t = ax.text(v, i, sizes2[i], color="black")
            bb = t.get_window_extent(renderer=r)
            bb = bb.transformed(transf)
            w_text = bb.width
            x_text = bb.x0
            y_text = bb.y0
            if v >= 0:
                t.set_position((x_text + w_text * 0.15, y_text))
            else:
                t.set_position((x_text - w_text, y_text))

    ax.grid(linestyle=":")

    if b is None:
        ax.legend(["Positive Weights", "Negative Weights"])
    else:
        ax.legend(["Positive Weights", "Negative Weights", b])

    if kind == "v":
        ax.axhline(y=0, xmin=0, xmax=1, color="gray", label=False)
    elif kind == "h":
        ax.axvline(x=0, ymin=0, ymax=1, color="gray", label=False)

    try:
        fig.set_layout_engine(layout="constrained")
    except:
        pass

    return ax


def plot_frontier_area(
    w_frontier, nrow=25, cmap="tab20", n_colors=20, height=6, width=10, ax=None
):
    r"""
    Create a chart with the asset composition of the efficient frontier.

    Parameters
    ----------
    w_frontier : DataFrame
        Weights of portfolios in the efficient frontier.
    nrow : int, optional
        Number of rows of the legend. The default is 25.
    cmap : cmap, optional
        Color scale used to plot each asset weight.
        The default is 'tab20'.
    n_colors : int, optional
        Number of distinct colors per color cycle. If there are more assets
        than n_colors, the chart is going to start to repeat the color cycle.
        The default is 20.
    height : float, optional
        Height of the image in inches. The default is 6.
    width : float, optional
        Width of the image in inches. The default is 10.
    ax : matplotlib axis, optional
        If provided, plot on this axis. The default is None.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    Returns
    -------
    ax :  matplotlib axis.
        Returns the Axes object with the plot for further tweaking.

    Example
    -------
    ::

        ax = rp.plot_frontier_area(w_frontier=ws,
                                   cmap="tab20",
                                   height=6,
                                   width=10,
                                   ax=None)

    .. image:: images/Area_Frontier.png


    """

    if not isinstance(w_frontier, pd.DataFrame):
        raise ValueError("w_frontier must be a DataFrame.")

    index = w_frontier.index.tolist()
    columns = w_frontier.columns.tolist()
    if not all(isinstance(x, str) for x in index):
        if not all(isinstance(x, str) for x in columns):
            raise ValueError("w_frontier index must be the names of assets.")
        else:
            w_frontier_ = w_frontier.T.copy()
    else:
        w_frontier_ = w_frontier.copy()

    columns = w_frontier_.columns.tolist()
    if not (
        all(isinstance(x, int) for x in columns)
        or all(isinstance(x, float) for x in columns)
    ):
        raise ValueError(
            "w_frontier columns must be the number of the point in the efficient frontier."
        )

    columns = list(range(len(columns)))

    if ax is None:
        fig = plt.gcf()
        ax = fig.gca()
        fig.set_figwidth(width)
        fig.set_figheight(height)
    else:
        fig = ax.get_figure()

    ax.set_title("Efficient Frontier's Assets Structure")

    colormap = cm.get_cmap(cmap)
    colormap = colormap(np.linspace(0, 1, int(n_colors)))

    if cmap == "gist_rainbow":
        colormap = colormap[::-1]

    cycle = plt.cycler("color", colormap)
    ax.set_prop_cycle(cycle)

    ax.stackplot(columns, w_frontier_, labels=index, alpha=0.7, edgecolor="black")

    ax.set_ylim(0, 1)
    ax.set_xlim(0, len(columns) - 1)

    ticks_loc = ax.get_yticks().tolist()
    ax.set_yticks(ax.get_yticks().tolist())
    ax.set_yticklabels(["{:3.2%}".format(x) for x in ticks_loc])
    ax.grid(linestyle=":")

    n = int(np.ceil(len(index) / nrow))

    ax.legend(index, loc="center left", bbox_to_anchor=(1, 0.5), ncol=n)

    try:
        fig.set_layout_engine(layout="constrained")
    except:
        pass

    return ax


def plot_risk_con(
    w,
    returns,
    cov=None,
    asset_classes=None,
    classes_col=None,
    rm="MV",
    rf=0,
    alpha=0.05,
    a_sim=100,
    beta=None,
    b_sim=None,
    kappa=0.30,
    kappa_g=None,
    solver="CLARABEL",
    percentage=False,
    erc_line=True,
    color="tab:blue",
    erc_linecolor="r",
    height=6,
    width=10,
    t_factor=252,
    ax=None,
):
    r"""
    Create a chart with the risk contribution per asset of the portfolio.

    Parameters
    ----------
    w : DataFrame or Series of shape (n_assets, 1)
        Portfolio weights, where n_assets is the number of assets.
    returns : DataFrame of shape (n_samples, n_assets)
        Assets returns DataFrame, where n_samples is the number of
        observations and n_assets is the number of assets.
    cov : DataFrame of shape (n_assets, n_assets)
        Covariance matrix, where n_assets is the number of assets.
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
    rm : str, optional
        Risk measure used to estimate risk contribution.
        The default is 'MV'. Possible values are:

        - 'MV': Standard Deviation.
        - 'KT': Square Root Kurtosis.
        - 'MAD': Mean Absolute Deviation.
        - 'GMD': Gini Mean Difference.
        - 'MSV': Semi Standard Deviation.
        - 'SKT': Square Root Semi Kurtosis.
        - 'FLPM': First Lower Partial Moment (Omega Ratio).
        - 'SLPM': Second Lower Partial Moment (Sortino Ratio).
        - 'CVaR': Conditional Value at Risk.
        - 'TG': Tail Gini.
        - 'EVaR': Entropic Value at Risk.
        - 'RLVaR': Relativistic Value at Risk.
        - 'WR': Worst Realization (Minimax).
        - 'CVRG': CVaR range of returns.
        - 'TGRG': Tail Gini range of returns.
        - 'EVRG': EVaR range of returns.
        - 'RVRG': RLVaR range of returns.
        - 'RG': Range of returns.
        - 'MDD': Maximum Drawdown of uncompounded cumulative returns (Calmar Ratio).
        - 'ADD': Average Drawdown of uncompounded cumulative returns.
        - 'CDaR': Conditional Drawdown at Risk of uncompounded cumulative returns.
        - 'EDaR': Entropic Drawdown at Risk of uncompounded cumulative returns.
        - 'RLDaR': Relativistic Drawdown at Risk of uncompounded cumulative returns.
        - 'UCI': Ulcer Index of uncompounded cumulative returns.

    rf : float, optional
        Risk free rate or minimum acceptable return. The default is 0.
    alpha : float, optional
        Significance level of VaR, CVaR, Tail Gini, EVaR, RLVaR, CDaR, EDaR and RLDaR. The default is 0.05.
    a_sim : float, optional
        Number of CVaRs used to approximate Tail Gini of losses. The default is 100.
    beta : float, optional
        Significance level of CVaR and Tail Gini of gains. If None it duplicates alpha value.
        The default is None.
    b_sim : float, optional
        Number of CVaRs used to approximate Tail Gini of gains. If None it duplicates a_sim value.
        The default is None.
    kappa : float, optional
        Deformation parameter of RLVaR and RLDaR for losses, must be between 0 and 1. The default is 0.30.
    kappa_g : float, optional
        Deformation parameter of RLVaR for gains, must be between 0 and 1.
        The default is None.
    solver: str, optional
        Solver available for CVXPY that supports power cone programming. Used to calculate RLVaR and RLDaR.
        The default value is 'CLARABEL'.
    percentage : bool, optional
        If risk contribution per asset is expressed as percentage or as a value. The default is False.
    erc_line : bool, optional
        If equal risk contribution line is plotted.
        The default is False.
    color : str, optional
        Color used to plot each asset risk contribution.
        The default is 'tab:blue'.
    erc_linecolor : str, optional
        Color used to plot equal risk contribution line.
        The default is 'r'.
    height : float, optional
        Height of the image in inches. The default is 6.
    width : float, optional
        Width of the image in inches. The default is 10.
    t_factor : float, optional
        Factor used to annualize expected return and expected risks for
        risk measures based on returns (not drawdowns). The default is 252.

        .. math::

            \begin{align}
            \text{Annualized Return} & = \text{Return} \, \times \, \text{t_factor} \\
            \text{Annualized Risk} & = \text{Risk} \, \times \, \sqrt{\text{t_factor}}
            \end{align}

    ax : matplotlib axis, optional
        If provided, plot on this axis. The default is None.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    Returns
    -------
    ax :  matplotlib axis.
        Returns the Axes object with the plot for further tweaking.

    Example
    -------
    ::

        ax = rp.plot_risk_con(w=w2,
                              cov=cov,
                              returns=Y,
                              rm=rm,
                              rf=0,
                              alpha=0.05,
                              color="tab:blue",
                              height=6,
                              width=10,
                              t_factor=252,
                              ax=None)

    .. image:: images/Risk_Con.png


    """

    if not isinstance(returns, pd.DataFrame):
        raise ValueError("returns must be a DataFrame.")

    if not isinstance(w, pd.DataFrame):
        if isinstance(w, pd.Series):
            w_ = w.to_frame()
        else:
            raise ValueError("w must be a one column DataFrame or Series")
    else:
        if w.shape[0] == 1:
            w_ = w.T.copy()
        elif w.shape[1] == 1:
            w_ = w.copy()
        else:
            raise ValueError("w must be a one column DataFrame or Series")

    if returns.columns.tolist() != w_.index.tolist():
        raise ValueError("returns and w must have same assets.")

    if cov is None:
        cov = returns.cov()
    elif isinstance(cov, pd.DataFrame):
        if cov.index.tolist() != cov.columns.tolist():
            raise ValueError(
                "cov must be a square DataFrame with samen labels in indexes and columns."
            )
        elif cov.index.tolist() != w_.index.tolist():
            raise ValueError("cov and w must have the same assets.")
    else:
        raise ValueError("cov must be a square DataFrame.")

    if asset_classes is not None and classes_col is not None:
        if not isinstance(asset_classes, pd.DataFrame):
            raise ValueError("asset_classes must be a DataFrame")
        else:
            if asset_classes.shape[1] < 2:
                raise ValueError("asset_classes must have at least two columns")

            classes = asset_classes.columns.tolist()

            if isinstance(classes_col, str) and classes_col in classes:
                A = asset_classes.loc[:, classes_col].to_frame()
                col = A.columns.to_list()[0]
            elif isinstance(classes_col, int) and classes[classes_col] in classes:
                A = asset_classes.iloc[:, classes_col].to_frame()
                col = A.columns.to_list()[0]
            else:
                raise ValueError(
                    "classes_col must be a valid column or column position of asset_classes"
                )

    if beta is None:
        beta = alpha
    if b_sim is None:
        b_sim = a_sim
    if kappa_g is None:
        kappa_g = kappa

    if ax is None:
        fig = plt.gcf()
        ax = fig.gca()
        fig.set_figwidth(width)
        fig.set_figheight(height)
    else:
        fig = ax.get_figure()

    item = rmeasures.index(rm)
    if rm in [
        "CVaR",
        "TG",
        "EVaR",
        "RLVaR",
        "CVRG",
        "TGRG",
        "EVRG",
        "RVRG",
        "CDaR",
        "EDaR",
        "RLDaR",
    ]:
        title = "Risk (" + rm_names[item] + " $\\alpha = $" + "{0:.2%}".format(alpha)
    else:
        title = "Risk (" + rm_names[item]
    if rm in ["CVRG", "TGRG", "EVRG", "RVRG"]:
        title += ", $\\beta = $" + "{0:.2%}".format(beta)
    if rm in ["RLVaR", "RLDaR", "RVRG"]:
        title += ", $\\kappa = $" + "{0:.2}".format(kappa)
    if rm in ["RLVaR", "RLDaR", "RVRG"]:
        title += ", $\\kappa = $" + "{0:.2}".format(kappa_g)

    title += ") Contribution per Asset"
    if percentage:
        title += " (%)"
    ax.set_title(r"{}".format(title))

    X = w_.index.tolist()

    RC = rk.Risk_Contribution(
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
    )

    if rm not in ["MDD", "ADD", "CDaR", "EDaR", "RLDaR", "UCI"]:
        RC = RC * t_factor**0.5

    if percentage:
        RC = RC / np.sum(RC)

    if asset_classes is not None and classes_col is not None:
        A = asset_classes.copy()
        B = pd.DataFrame(RC, index=X)
        A = pd.merge(A, B, left_on=classes[0], right_index=True, how="left")
        A = A.groupby([col]).sum()[0]
        X = A.index.tolist()
        RC = A.to_numpy()

    ax.bar(X, RC, alpha=0.7, color=color, edgecolor="black")

    ax.set_xlim(-0.5, len(X) - 0.5)
    ax.tick_params(axis="x", labelrotation=90)

    ticks_loc = ax.get_yticks().tolist()
    ax.set_yticks(ax.get_yticks())
    ax.set_yticklabels(["{:3.4%}".format(x) for x in ticks_loc])
    ax.grid(linestyle=":")

    if erc_line:
        if percentage:
            erc = 1 / len(RC)
        else:
            erc = rk.Sharpe_Risk(
                returns=returns,
                w=w_,
                cov=cov,
                rm=rm,
                rf=rf,
                alpha=alpha,
                a_sim=a_sim,
                beta=beta,
                b_sim=b_sim,
                kappa=kappa,
                solver=solver,
            )

            if rm not in ["MDD", "ADD", "CDaR", "EDaR", "RLDaR", "UCI"]:
                erc = erc / len(RC) * t_factor**0.5
            else:
                erc = erc / len(RC)

        ax.axhline(y=erc, color=erc_linecolor, linestyle="-")

    try:
        fig.set_layout_engine(layout="constrained")
    except:
        pass

    return ax


def plot_factor_risk_con(
    w,
    returns,
    factors,
    cov=None,
    B=None,
    const=True,
    rm="MV",
    rf=0,
    alpha=0.05,
    a_sim=100,
    beta=None,
    b_sim=None,
    kappa=0.30,
    kappa_g=None,
    solver="CLARABEL",
    feature_selection="stepwise",
    stepwise="Forward",
    criterion="pvalue",
    threshold=0.05,
    n_components=0.95,
    percentage=False,
    erc_line=True,
    color="tab:orange",
    erc_linecolor="r",
    height=6,
    width=10,
    t_factor=252,
    ax=None,
):
    r"""
    Create a chart with the risk contribution per risk factor of the portfolio.

    Parameters
    ----------
    w : DataFrame  or Series of shape (n_assets, 1)
        Portfolio weights, where n_assets is the number of assets.
    cov : DataFrame of shape (n_assets, n_assets)
        Covariance matrix, where n_assets is the number of assets.
    returns : DataFrame of shape (n_samples, n_assets)
        Assets returns DataFrame, where n_samples is the number of
        observations and n_assets is the number of assets.
    factors : DataFrame or nd-array of shape (n_samples, n_factors)
        Risk factors returns DataFrame, where n_samples is the number of samples and
        n_factors is the number of factors.
    B : DataFrame of shape (n_assets, n_factors), optional
        Loadings matrix, where n_assets is the number assets and n_factors is
        the number of risk factors. If is not specified, is estimated using
        stepwise regression. The default is None.
    const : bool, optional
        Indicate if the loadings matrix has a constant.
        The default is False.
    rm : str, optional
        Risk measure used to estimate risk contribution.
        The default is 'MV'. Possible values are:

        - 'MV': Standard Deviation.
        - 'KT': Square Root Kurtosis.
        - 'MAD': Mean Absolute Deviation.
        - 'GMD': Gini Mean Difference.
        - 'MSV': Semi Standard Deviation.
        - 'SKT': Square Root Semi Kurtosis.
        - 'FLPM': First Lower Partial Moment (Omega Ratio).
        - 'SLPM': Second Lower Partial Moment (Sortino Ratio).
        - 'CVaR': Conditional Value at Risk.
        - 'TG': Tail Gini.
        - 'EVaR': Entropic Value at Risk.
        - 'RLVaR': Relativistic Value at Risk.
        - 'WR': Worst Realization (Minimax).
        - 'CVRG': CVaR range of returns.
        - 'TGRG': Tail Gini range of returns.
        - 'EVRG': EVaR range of returns.
        - 'RVRG': RLVaR range of returns.
        - 'RG': Range of returns.
        - 'MDD': Maximum Drawdown of uncompounded cumulative returns (Calmar Ratio).
        - 'ADD': Average Drawdown of uncompounded cumulative returns.
        - 'CDaR': Conditional Drawdown at Risk of uncompounded cumulative returns.
        - 'EDaR': Entropic Drawdown at Risk of uncompounded cumulative returns.
        - 'RLDaR': Relativistic Drawdown at Risk of uncompounded cumulative returns.
        - 'UCI': Ulcer Index of uncompounded cumulative returns.

    rf : float, optional
        Risk free rate or minimum acceptable return. The default is 0.
    alpha : float, optional
        Significance level of VaR, CVaR, Tail Gini, EVaR, RLVaR, CDaR, EDaR and RLDaR. The default is 0.05.
    a_sim : float, optional
        Number of CVaRs used to approximate Tail Gini of losses. The default is 100.
    beta : float, optional
        Significance level of CVaR and Tail Gini of gains. If None it duplicates alpha value.
        The default is None.
    b_sim : float, optional
        Number of CVaRs used to approximate Tail Gini of gains. If None it duplicates a_sim value.
        The default is None.
    kappa : float, optional
        Deformation parameter of RLVaR and RLDaR for losses, must be between 0 and 1. The default is 0.30.
    kappa_g : float, optional
        Deformation parameter of RLVaR for gains, must be between 0 and 1.
        The default is None.
    solver: str, optional
        Solver available for CVXPY that supports power cone programming. Used to calculate RLVaR and RLDaR.
        The default value is 'CLARABEL'.
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
    percentage : bool, optional
        If risk contribution per asset is expressed as percentage or as a value. The default is False.
    erc_line : bool, optional
        If equal risk contribution line is plotted.
        The default is False.
    color : str, optional
        Color used to plot each asset risk contribution.
        The default is 'tab:orange'.
    erc_linecolor : str, optional
        Color used to plot equal risk contribution line.
        The default is 'r'.
    height : float, optional
        Height of the image in inches. The default is 6.
    width : float, optional
        Width of the image in inches. The default is 10.
    t_factor : float, optional
        Factor used to annualize expected return and expected risks for
        risk measures based on returns (not drawdowns). The default is 252.

        .. math::

            \begin{align}
            \text{Annualized Return} & = \text{Return} \, \times \, \text{t_factor} \\
            \text{Annualized Risk} & = \text{Risk} \, \times \, \sqrt{\text{t_factor}}
            \end{align}

    ax : matplotlib axis, optional
        If provided, plot on this axis. The default is None.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    Returns
    -------
    ax :  matplotlib axis.
        Returns the Axes object with the plot for further tweaking.

    Example
    -------
    ::

        ax = rp.plot_factor_risk_con(w=w3,
                                     cov=cov,
                                     returns=Y,
                                     factors=X,
                                     B=None,
                                     const=True,
                                     rm=rm,
                                     rf=0,
                                     feature_selection="stepwise",
                                     stepwise="Forward",
                                     criterion="pvalue",
                                     threshold=0.05,
                                     height=6,
                                     width=10,
                                     t_factor=252,
                                     ax=None)

    .. image:: images/Risk_Con_RF.png

    ::

        ax = rp.plot_factor_risk_con(w=w4,
                                     cov=cov,
                                     returns=Y,
                                     factors=X,
                                     B=None,
                                     const=True,
                                     rm=rm,
                                     rf=0,
                                     feature_selection="PCR",
                                     n_components=0.95,
                                     height=6,
                                     width=10,
                                     t_factor=252,
                                     ax=None)

    .. image:: images/Risk_Con_PC.png

    """

    if not isinstance(w, pd.DataFrame):
        if isinstance(w, pd.Series):
            w_ = w.to_frame()
        else:
            raise ValueError("w must be a one column DataFrame or Series")
    else:
        if w.shape[0] == 1:
            w_ = w.T.copy()
        elif w.shape[1] == 1:
            w_ = w.copy()
        else:
            raise ValueError("w must be a one column DataFrame or Series")

    if returns.columns.tolist() != w_.index.tolist():
        raise ValueError("returns and w must have same assets.")

    if cov is None:
        cov = returns.cov()
    elif isinstance(cov, pd.DataFrame):
        if cov.index.tolist() != cov.columns.tolist():
            raise ValueError(
                "cov must be a square DataFrame with samen labels in indexes and columns."
            )
        elif cov.index.tolist() != w_.index.tolist():
            raise ValueError("cov and w must have the same assets.")
    else:
        raise ValueError("cov must be a square DataFrame.")

    if beta is None:
        beta = alpha
    if b_sim is None:
        b_sim = a_sim
    if kappa_g is None:
        kappa_g = kappa

    if ax is None:
        fig = plt.gcf()
        ax = fig.gca()
        fig.set_figwidth(width)
        fig.set_figheight(height)
    else:
        fig = ax.get_figure()

    item = rmeasures.index(rm)
    if rm in [
        "CVaR",
        "TG",
        "EVaR",
        "RLVaR",
        "CVRG",
        "TGRG",
        "EVRG",
        "RVRG",
        "CDaR",
        "EDaR",
        "RLDaR",
    ]:
        title = "Risk (" + rm_names[item] + " $\\alpha = $" + "{0:.2%}".format(alpha)
    else:
        title = "Risk (" + rm_names[item]
    if rm in ["CVRG", "TGRG", "EVRG", "RVRG"]:
        title += ", $\\beta = $" + "{0:.2%}".format(beta)
    if rm in ["RLVaR", "RLDaR", "RVRG"]:
        title += ", $\\kappa = $" + "{0:.2}".format(kappa)
    if rm in ["RLVaR", "RLDaR", "RVRG"]:
        title += ", $\\kappa = $" + "{0:.2}".format(kappa_g)

    title += ") Contribution per "
    if feature_selection == "PCR":
        title += "Principal Component (PC)"
    else:
        title += "Risk Factor"
    if percentage:
        title += " (%)"
    ax.set_title(r"{}".format(title))

    RC_F = rk.Factors_Risk_Contribution(
        w=w_,
        returns=returns,
        factors=factors,
        cov=cov,
        B=B,
        const=const,
        rm=rm,
        rf=rf,
        alpha=alpha,
        a_sim=a_sim,
        beta=beta,
        b_sim=b_sim,
        kappa=kappa,
        solver=solver,
        feature_selection=feature_selection,
        stepwise=stepwise,
        criterion=criterion,
        threshold=threshold,
        n_components=n_components,
    )

    if feature_selection == "PCR":
        n = len(RC_F) - 1
        X = ["PC " + str(i) for i in range(1, n + 1)]
    else:
        X = factors.columns.tolist()
    X.append("Others")

    if rm not in ["MDD", "ADD", "CDaR", "EDaR", "RLDaR", "UCI"]:
        RC_F = RC_F * t_factor**0.5

    if percentage:
        RC_F = RC_F / np.sum(RC_F)

    ax.bar(X, RC_F, alpha=0.7, color=color, edgecolor="black")

    ax.set_xlim(-0.5, len(X) - 0.5)
    ax.tick_params(axis="x", labelrotation=90)

    ticks_loc = ax.get_yticks().tolist()
    ax.set_yticks(ax.get_yticks())
    ax.set_yticklabels(["{:3.4%}".format(x) for x in ticks_loc])
    ax.grid(linestyle=":")

    if erc_line:
        if percentage:
            erc = 1 / len(RC_F)
        else:
            erc = rk.Sharpe_Risk(
                returns=returns,
                w=w_,
                cov=cov,
                rm=rm,
                rf=rf,
                alpha=alpha,
                a_sim=a_sim,
                beta=beta,
                b_sim=b_sim,
                kappa=kappa,
                solver=solver,
            )

            if rm not in ["MDD", "ADD", "CDaR", "EDaR", "RLDaR", "UCI"]:
                erc = erc / len(RC_F) * t_factor**0.5
            else:
                erc = erc / len(RC_F)

        ax.axhline(y=erc, color=erc_linecolor, linestyle="-")

    try:
        fig.set_layout_engine(layout="constrained")
    except:
        pass

    return ax


def plot_hist(
    returns,
    w,
    alpha=0.05,
    a_sim=100,
    kappa=0.30,
    solver="CLARABEL",
    bins=50,
    height=6,
    width=10,
    ax=None,
):
    r"""
    Create a histogram of portfolio returns with the risk measures.

    Parameters
    ----------
    returns : DataFrame of shape (n_samples, n_assets)
        Assets returns DataFrame, where n_samples is the number of
        observations and n_assets is the number of assets.
    w : DataFrame or Series of shape (n_assets, 1)
        Portfolio weights, where n_assets is the number of assets.
    alpha : float, optional
        Significance level of VaR, CVaR, EVaR, RLVaR and Tail Gini. The default is 0.05.
    a_sim : float, optional
        Number of CVaRs used to approximate Tail Gini of losses. The default is 100.
    kappa : float, optional
        Deformation parameter of RLVaR and RLDaR, must be between 0 and 1. The default is 0.30.
    solver: str, optional
        Solver available for CVXPY that supports power cone programming. Used to calculate RLVaR and RLDaR.
        The default value is 'CLARABEL'.
    bins : float, optional
        Number of bins of the histogram. The default is 50.
    height : float, optional
        Height of the image in inches. The default is 6.
    width : float, optional
        Width of the image in inches. The default is 10.
    ax : matplotlib axis, optional
        If provided, plot on this axis. The default is None.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    Returns
    -------
    ax : matplotlib axis.
        Returns the Axes object with the plot for further tweaking.

    Example
    -------
    ::

        ax = rp.plot_hist(returns=Y,
                          w=w1,
                          alpha=0.05,
                          bins=50,
                          height=6,
                          width=10,
                          ax=None)

    .. image:: images/Histogram.png


    """

    if not isinstance(returns, pd.DataFrame):
        raise ValueError("returns must be a DataFrame")

    if not isinstance(w, pd.DataFrame):
        if isinstance(w, pd.Series):
            w_ = w.to_frame()
        else:
            raise ValueError("w must be a one column DataFrame or Series")
    else:
        if w.shape[0] == 1:
            w_ = w.T.copy()
        elif w.shape[1] == 1:
            w_ = w.copy()
        else:
            raise ValueError("w must be a one column DataFrame or Series")

    if returns.columns.tolist() != w_.index.tolist():
        raise ValueError("returns and w must have same assets.")

    if ax is None:
        fig = plt.gcf()
        ax = fig.gca()
        fig.set_figwidth(width)
        fig.set_figheight(height)
    else:
        fig = ax.get_figure()

    a = returns.to_numpy() @ w_.to_numpy()
    ax.set_title("Portfolio Returns Histogram")
    n, bins1, patches = ax.hist(
        a, int(bins), density=1, edgecolor="skyblue", color="skyblue", alpha=0.5
    )
    mu = np.mean(a)
    sigma = np.std(a, axis=0, ddof=1).item()
    risk = [
        mu,
        mu - sigma,
        mu - rk.MAD(a),
        mu - rk.GMD(a),
        -rk.VaR_Hist(a, alpha),
        -rk.CVaR_Hist(a, alpha),
        -rk.TG(a, alpha, a_sim),
        -rk.EVaR_Hist(a, alpha)[0],
        -rk.RLVaR_Hist(a, alpha, kappa, solver),
        -rk.WR(a),
    ]
    label = [
        "Mean: " + "{0:.2%}".format(risk[0]),
        "Mean - Std. Dev.("
        + "{0:.2%}".format(-risk[1] + mu)
        + "): "
        + "{0:.2%}".format(risk[1]),
        "Mean - MAD("
        + "{0:.2%}".format(-risk[2] + mu)
        + "): "
        + "{0:.2%}".format(risk[2]),
        "Mean - GMD("
        + "{0:.2%}".format(-risk[3] + mu)
        + "): "
        + "{0:.2%}".format(risk[3]),
        "{0:.2%}".format((1 - alpha)) + " Confidence VaR: " + "{0:.2%}".format(risk[4]),
        "{0:.2%}".format((1 - alpha))
        + " Confidence CVaR: "
        + "{0:.2%}".format(risk[5]),
        "{0:.2%}".format((1 - alpha))
        + " Confidence Tail Gini: "
        + "{0:.2%}".format(risk[6]),
        "{0:.2%}".format((1 - alpha))
        + " Confidence EVaR: "
        + "{0:.2%}".format(risk[7]),
        "{0:.1%}".format((1 - alpha))
        + " Confidence RLVaR("
        + "{0:.3}".format(kappa)
        + "): "
        + "{0:.2%}".format(risk[8]),
        "Worst Realization: " + "{0:.2%}".format(risk[9]),
    ]
    color = [
        "b",
        "r",
        "fuchsia",
        "navy",
        "darkorange",
        "limegreen",
        "mediumvioletred",
        "dodgerblue",
        "slateblue",
        "darkgrey",
    ]

    for i, j, k in zip(risk, label, color):
        ax.axvline(x=i, color=k, linestyle="-", label=j)

    # add a 'best fit' line
    y = (1 / (np.sqrt(2 * np.pi) * sigma)) * np.exp(
        -0.5 * (1 / sigma * (bins1 - mu)) ** 2
    )
    ax.plot(
        bins1,
        y,
        "--",
        color="orange",
        label=r"Normal: $\mu=$"
        + r"{0:.2%},".format(mu)
        + r" $\sigma=$"
        + r"{0:.2%}".format(sigma),
    )

    factor = (np.max(a) - np.min(a)) / bins

    ax.xaxis.set_major_locator(plt.AutoLocator())
    ticks_loc = ax.get_xticks().tolist()
    ticks_loc = [round(i, 8) for i in ticks_loc]
    ticks_loc = ticks_loc + [-i for i in ticks_loc[::-1]]
    ticks_loc = list(set(ticks_loc))
    ticks_loc.sort()
    ax.set_xticks(np.array(ticks_loc))
    ax.set_xticklabels(["{:3.2%}".format(x) for x in ticks_loc])
    ticks_loc = ax.get_yticks().tolist()
    ax.set_yticks(ax.get_yticks())
    ax.set_yticklabels(["{:3.2%}".format(x * factor) for x in ticks_loc])
    ax.legend(loc="upper right")  # , fontsize = 'x-small')
    ax.grid(linestyle=":")
    ax.set_ylabel("Probability Density")

    try:
        fig.set_layout_engine(layout="constrained")
    except:
        pass

    return ax


def plot_range(
    returns,
    w,
    alpha=0.05,
    a_sim=100,
    beta=None,
    b_sim=None,
    kappa=0.3,
    kappa_g=0.3,
    solver="CLARABEL",
    bins=50,
    height=6,
    width=10,
    ax=None,
):
    r"""
    Create a histogram of portfolio returns with the range risk measures.

    Parameters
    ----------
    returns : DataFrame of shape (n_samples, n_assets)
        Assets returns DataFrame, where n_samples is the number of
        observations and n_assets is the number of assets.
    w : DataFrame or Series of shape (n_assets, 1)
        Portfolio weights, where n_assets is the number of assets.
    alpha : float, optional
        Significance level of CVaR and Tail Gini of losses.
        The default is 0.05.
    a_sim : float, optional
        Number of CVaRs used to approximate Tail Gini of losses. The default is 100.
    beta : float, optional
        Significance level of CVaR and Tail Gini of gains. If None it duplicates alpha value.
        The default is None.
    b_sim : float, optional
        Number of CVaRs used to approximate Tail Gini of gains. If None it duplicates a_sim value.
        The default is None.
    bins : float, optional
        Number of bins of the histogram. The default is 50.
    height : float, optional
        Height of the image in inches. The default is 6.
    width : float, optional
        Width of the image in inches. The default is 10.
    ax : matplotlib axis, optional
        If provided, plot on this axis. The default is None.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    Returns
    -------
    ax : matplotlib axis.
        Returns the Axes object with the plot for further tweaking.

    Example
    -------
    ::

        ax = rp.plot_range(returns=Y,
                           w=w1,
                           alpha=0.05,
                           a_sim=100,
                           beta=None,
                           b_sim=None,
                           bins=50,
                           height=6,
                           width=10,
                           ax=None)

    .. image:: images/Range.png


    """

    if not isinstance(returns, pd.DataFrame):
        raise ValueError("returns must be a DataFrame")

    if not isinstance(w, pd.DataFrame):
        if isinstance(w, pd.Series):
            w_ = w.to_frame()
        else:
            raise ValueError("w must be a one column DataFrame or Series")
    else:
        if w.shape[0] == 1:
            w_ = w.T.copy()
        elif w.shape[1] == 1:
            w_ = w.copy()
        else:
            raise ValueError("w must be a one column DataFrame or Series")

    if returns.columns.tolist() != w_.index.tolist():
        raise ValueError("returns and w must have same assets.")

    if beta is None:
        beta = alpha
    if b_sim is None:
        b_sim = a_sim
    if kappa_g is None:
        kappa_g = kappa

    if ax is None:
        fig = plt.gcf()
        ax = fig.gca()
        fig.set_figwidth(width)
        fig.set_figheight(height)
    else:
        fig = ax.get_figure()

    a = returns.to_numpy() @ w_.to_numpy()
    ax.set_title("Portfolio Returns Range Risk Measures")

    df = dict(
        risk=[
            "VaR Range",
            "CVaR Range",
            "Tail Gini Range",
            "EVaR Range",
            "RLVaR Range",
            "Range",
        ],
        lower=[],
        upper=[],
    )

    df["lower"].append(-rk.VaR_Hist(a, alpha=alpha))
    df["lower"].append(-rk.CVaR_Hist(a, alpha=alpha))
    df["lower"].append(-rk.TG(a, alpha=alpha, a_sim=a_sim))
    df["lower"].append(-rk.EVaR_Hist(a, alpha=alpha, solver=solver)[0])
    df["lower"].append(-rk.RLVaR_Hist(a, alpha=alpha, kappa=kappa, solver=solver))
    df["lower"].append(np.min(a))
    df["upper"].append(rk.VaR_Hist(-a, alpha=beta))
    df["upper"].append(rk.CVaR_Hist(-a, alpha=beta))
    df["upper"].append(rk.TG(-a, alpha=beta, a_sim=b_sim))
    df["upper"].append(rk.EVaR_Hist(-a, alpha=beta, solver=solver)[0])
    df["upper"].append(
        rk.RLVaR_Hist(-a * 100, alpha=beta, kappa=kappa_g, solver=solver) / 100
    )
    df["upper"].append(-np.min(-a))

    df = pd.DataFrame(df)
    df.set_index("risk", inplace=True)

    # Func to draw line segment
    def newline(p1, p2, color="black"):
        ax = fig.gca()
        l = mlines.Line2D([p1[0], p2[0]], [p1[1], p2[1]], color=color)
        ax.add_line(l)
        return l

    n, _, _ = ax.hist(a, bins=int(bins), density=True, color="darkgrey", alpha=0.3)

    risk = [
        rk.VRG(a, alpha=alpha, beta=beta),
        rk.CVRG(a, alpha=alpha, beta=beta),
        rk.TGRG(a, alpha=alpha, a_sim=a_sim, beta=beta, b_sim=b_sim),
        rk.EVRG(a, alpha=alpha, beta=beta, solver=solver),
        rk.RVRG(a, alpha=alpha, beta=beta, kappa=kappa, kappa_g=kappa_g, solver=solver),
        rk.RG(a),
    ]

    label = [
        "VaR Range ("
        + "{0:.1%}".format((1 - alpha))
        + ", "
        + "{0:.1%}".format((1 - beta))
        + "): "
        + "{0:.2%}".format(risk[0]),
        "CVaR Range ("
        + "{0:.1%}".format((1 - alpha))
        + ", "
        + "{0:.1%}".format((1 - beta))
        + "): "
        + "{0:.2%}".format(risk[1]),
        "Tail Gini Range ("
        + "{0:.1%}".format((1 - alpha))
        + ", "
        + "{0:.1%}".format((1 - beta))
        + "): "
        + "{0:.2%}".format(risk[2]),
        "EVaR Range ("
        + "{0:.1%}".format((1 - alpha))
        + ", "
        + "{0:.1%}".format((1 - beta))
        + "): "
        + "{0:.2%}".format(risk[3]),
        "RlVaR Range ("
        + "{0:.1%}".format((1 - alpha))
        + ", "
        + "{0:.1%}".format((1 - beta))
        + ", "
        + "{0:.1%}".format(kappa)
        + ", "
        + "{0:.1%}".format(kappa_g)
        + "): "
        + "{0:.2%}".format(risk[4]),
        "Range :" + "{0:.2%}".format(risk[5]),
    ]

    colors = [
        "darkorange",
        "limegreen",
        "mediumvioletred",
        "dodgerblue",
        "slateblue",
        "fuchsia",
    ]

    y_max = np.ceil(n.max())

    j = 1
    for i in df.index:
        x1 = df.loc[i, "lower"]
        x2 = df.loc[i, "upper"]
        y1 = (len(colors) + 1 - j) * y_max / 9
        ax.vlines(
            x=x1,
            ymin=0,
            ymax=y1,
            color=colors[j - 1],
            alpha=1,
            linewidth=1,
            linestyles="dashed",
        )
        ax.vlines(
            x=x2,
            ymin=0,
            ymax=y1,
            color=colors[j - 1],
            alpha=1,
            linewidth=1,
            linestyles="dashed",
        )
        ax.scatter(y=y1, x=x1, s=50, color=colors[j - 1], alpha=1, label=label[j - 1])
        ax.scatter(y=y1, x=x2, s=50, color=colors[j - 1], alpha=1)
        newline([x1, y1], [x2, y1], color=colors[j - 1])
        j += 1

    ax.set(ylim=(0, y_max))

    factor = (np.max(a) - np.min(a)) / bins
    ax.xaxis.set_major_locator(plt.AutoLocator())
    ticks_loc = ax.get_xticks().tolist()
    ax.set_xticks(ax.get_xticks())
    ax.set_xticklabels(["{:3.2%}".format(x) for x in ticks_loc])
    ticks_loc = ax.get_yticks().tolist()
    ax.set_yticks(ax.get_yticks())
    ax.set_yticklabels(["{:3.2%}".format(x * factor) for x in ticks_loc])
    ax.legend(loc="upper right")  # , fontsize = 'x-small')
    ax.grid(linestyle=":")
    ax.set_ylabel("Probability Density")

    try:
        fig.set_layout_engine(layout="constrained")
    except:
        pass

    return ax


def plot_drawdown(
    returns,
    w,
    alpha=0.05,
    kappa=0.30,
    solver="CLARABEL",
    height=8,
    width=10,
    height_ratios=[2, 3],
    ax=None,
):
    r"""
    Create a chart with the evolution of portfolio prices and drawdown.

    Parameters
    ----------
    returns : DataFrame of shape (n_samples, n_assets)
        Assets returns DataFrame, where n_samples is the number of
        observations and n_assets is the number of assets.
    w : DataFrame or Series of shape (n_assets, 1)
        Portfolio weights, where n_assets is the number of assets.
    alpha : float, optional
        Significance level of DaR, CDaR, EDaR and RLDaR. The default is 0.05.
    kappa : float, optional
        Deformation parameter of RLVaR and RLDaR, must be between 0 and 1. The default is 0.30.
    solver: str, optional
        Solver available for CVXPY that supports power cone programming. Used to calculate RLVaR and RLDaR.
        The default value is 'CLARABEL'.
    height : float, optional
        Height of the image in inches. The default is 8.
    width : float, optional
        Width of the image in inches. The default is 10.
    height_ratios : list or ndarray
        Defines the relative heights of the rows. Each row gets a relative
        height of height_ratios[i] / sum(height_ratios). The default value is
        [2,3].
    ax : matplotlib axis of size (2,1), optional
        If provided, plot on this axis. The default is None.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    Returns
    -------
    ax : np.array
        Returns the a np.array with Axes objects with plots for further tweaking.

    Example
    -------
    ::

        ax = rp.plot_drawdown(returns=Y,
                              w=w1,
                              alpha=0.05,
                              height=8,
                              width=10,
                              ax=None)

    .. image:: images/Drawdown.png


    """

    if not isinstance(returns, pd.DataFrame):
        raise ValueError("returns must be a DataFrame")

    if not isinstance(w, pd.DataFrame):
        if isinstance(w, pd.Series):
            w_ = w.to_frame()
        else:
            raise ValueError("w must be a one column DataFrame or Series")
    else:
        if w.shape[0] == 1:
            w_ = w.T.copy()
        elif w.shape[1] == 1:
            w_ = w.copy()
        else:
            raise ValueError("w must be a one column DataFrame or Series")

    if returns.columns.tolist() != w_.index.tolist():
        raise ValueError("returns and w must have same assets.")

    if ax is None:
        fig = plt.gcf()
        ax = fig.gca()
        ax.axis("off")
        gs = GridSpec(nrows=2, ncols=1, figure=fig, height_ratios=height_ratios)
        axes = []
        axes.append(fig.add_subplot(gs[0]))
        axes.append(fig.add_subplot(gs[1]))
        fig.set_figwidth(width)
        fig.set_figheight(height)
    else:
        if isinstance(ax, plt.Axes):
            ax.axis("off")
            fig = ax.get_figure()
            if hasattr(ax, "get_subplotspec"):
                subplot_spec = ax.get_subplotspec()
                gs0 = subplot_spec.get_gridspec()
                gs = GridSpecFromSubplotSpec(
                    nrows=2, ncols=1, height_ratios=height_ratios, subplot_spec=gs0[0]
                )
            else:
                gs = GridSpec(nrows=2, ncols=1, figure=fig, height_ratios=height_ratios)
            axes = []
            axes.append(fig.add_subplot(gs[0]))
            axes.append(fig.add_subplot(gs[1]))
        else:
            raise TypeError("ax must be a matplotlib axes object.")

    index = returns.index.tolist()
    a = returns.to_numpy() @ w_.to_numpy()
    prices = 1 + np.insert(a, 0, 0, axis=0)
    prices = np.cumprod(prices, axis=0)
    prices = np.ravel(prices).tolist()
    prices2 = 1 + np.array(np.cumsum(a, axis=0))
    prices2 = np.ravel(prices2).tolist()
    del prices[0]

    DD = []
    peak = -99999
    for i in range(0, len(prices)):
        if prices2[i] > peak:
            peak = prices2[i]
        DD.append((peak - prices2[i]))
    DD = -np.array(DD)
    titles = [
        "Historical Compounded Cumulative Returns",
        "Historical Uncompounded Drawdown",
    ]
    data = [prices, DD]
    color1 = ["b", "orange"]
    risk = [
        -rk.UCI_Abs(a),
        -rk.ADD_Abs(a),
        -rk.DaR_Abs(a, alpha),
        -rk.CDaR_Abs(a, alpha),
        -rk.EDaR_Abs(a, alpha)[0],
        -rk.RLDaR_Abs(a, alpha, kappa, solver),
        -rk.MDD_Abs(a),
    ]
    label = [
        "Ulcer Index: " + "{0:.2%}".format(risk[0]),
        "Average Drawdown: " + "{0:.2%}".format(risk[1]),
        "{0:.2%}".format((1 - alpha)) + " Confidence DaR: " + "{0:.2%}".format(risk[2]),
        "{0:.2%}".format((1 - alpha))
        + " Confidence CDaR: "
        + "{0:.2%}".format(risk[3]),
        "{0:.2%}".format((1 - alpha))
        + " Confidence EDaR: "
        + "{0:.2%}".format(risk[4]),
        "{0:.2%}".format((1 - alpha))
        + " Confidence RLDaR ("
        + "{0:.3}".format(kappa)
        + "): "
        + "{0:.2%}".format(risk[5]),
        "Maximum Drawdown: " + "{0:.2%}".format(risk[6]),
    ]
    color2 = ["b", "r", "fuchsia", "limegreen", "dodgerblue", "slateblue", "darkgrey"]

    j = 0

    ymin = np.min(DD) * 1.5

    locator = mdates.AutoDateLocator(minticks=5, maxticks=10)
    formatter = mdates.DateFormatter("%Y-%m")
    for i in axes:
        i.clear()
        i.plot_date(index, data[j], "-", color=color1[j])
        if j == 1:
            i.fill_between(index, 0, data[j], facecolor=color1[j], alpha=0.3)
            for k in range(0, len(risk)):
                i.axhline(y=risk[k], color=color2[k], linestyle="-", label=label[k])
            i.set_ylim(ymin, 0)
            i.legend(loc="lower right")  # , fontsize = 'x-small')
        i.set_title(titles[j])
        i.xaxis.set_major_locator(locator)
        i.xaxis.set_major_formatter(formatter)
        ticks_loc = i.get_yticks().tolist()
        i.set_yticks(i.get_yticks())
        i.set_yticklabels(["{:3.2%}".format(x) for x in ticks_loc])
        i.grid(linestyle=":")
        j = j + 1

    try:
        fig.set_layout_engine(layout="constrained")
    except:
        pass

    return ax


def plot_table(
    returns,
    w,
    MAR=0,
    alpha=0.05,
    a_sim=100,
    kappa=0.30,
    solver="CLARABEL",
    height=9,
    width=12,
    t_factor=252,
    ini_days=1,
    days_per_year=252,
    ax=None,
):
    r"""
    Create a table with information about risk measures and risk adjusted
    return ratios.

    Parameters
    ----------
    returns : DataFrame of shape (n_samples, n_assets)
        Assets returns DataFrame, where n_samples is the number of
        observations and n_assets is the number of assets.
    w : DataFrame or Series of shape (n_assets, 1)
        Portfolio weights, where n_assets is the number of assets.
    MAR: float, optional
        Minimum acceptable return.
    alpha : float, optional
        Significance level of VaR, CVaR, Tail Gini, EVaR, RLVaR, CDaR, EDaR and RLDaR.
        The default is 0.05.
    a_sim : float, optional
        Number of CVaRs used to approximate Tail Gini of losses. The default is 100.
    kappa : float, optional
        Deformation parameter of RLVaR and RLDaR, must be between 0 and 1. The default is 0.30.
    solver: str, optional
        Solver available for CVXPY that supports power cone programming. Used to calculate RLVaR and RLDaR.
        The default value is 'CLARABEL'.
    height : float, optional
        Height of the image in inches. The default is 9.
    width : float, optional
        Width of the image in inches. The default is 12.
    t_factor : float, optional
        Factor used to annualize expected return and expected risks for
        risk measures based on returns (not drawdowns). The default is 252.

        .. math::

            \begin{align}
            \text{Annualized Return} & = \text{Return} \, \times \, \text{t_factor} \\
            \text{Annualized Risk} & = \text{Risk} \, \times \, \sqrt{\text{t_factor}}
            \end{align}

    ini_days : float, optional
        If provided, it is the number of days of compounding for first return.
        It is used to calculate Compound Annual Growth Rate (CAGR). This value
        depend on assumptions used in t_factor, for example if data is monthly
        you can use 21 (252 days per year) or 30 (360 days per year). The
        default is 1 for daily returns.
    days_per_year: float, optional
        Days per year assumption. It is used to calculate Compound Annual
        Growth Rate (CAGR). Default value is 252 trading days per year.
    ax : matplotlib axis, optional
        If provided, plot on this axis. The default is None.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    Returns
    -------
    ax : matplotlib axis
        Returns the Axes object with the plot for further tweaking.

    Example
    -------
    ::

        ax = rp.plot_table(returns=Y,
                           w=w1,
                           MAR=0,
                           alpha=0.05,
                           ax=None)

    .. image:: images/Port_Table.png


    """
    if not isinstance(returns, pd.DataFrame):
        raise ValueError("returns must be a DataFrame")

    if not isinstance(w, pd.DataFrame):
        if isinstance(w, pd.Series):
            w_ = w.to_frame()
        else:
            raise ValueError("w must be a one column DataFrame or Series")
    else:
        if w.shape[0] == 1:
            w_ = w.T.copy()
        elif w.shape[1] == 1:
            w_ = w.copy()
        else:
            raise ValueError("w must be a one column DataFrame or Series")

    if returns.columns.tolist() != w_.index.tolist():
        raise ValueError("returns and w must have same assets.")

    if ax is None:
        fig = plt.gcf()
        ax = fig.gca()
        fig.set_figwidth(width)
        fig.set_figheight(height)
    else:
        fig = ax.get_figure()

    mu = returns.mean()
    cov = returns.cov()
    days = (returns.index[-1] - returns.index[0]).days + ini_days

    X = (returns.to_numpy() @ w_.to_numpy()).flatten()

    rowLabels = [
        "Profitability and Other Inputs",
        "Mean Return (1)",
        "Compound Annual Growth Rate (CAGR)",
        "Minimum Acceptable Return (MAR) (1)",
        "Significance Level",
        "",
        "Risk Measures based on Returns",
        "Standard Deviation (2)",
        "Mean Absolute Deviation (MAD) (2)",
        "Semi Standard Deviation (2)",
        "First Lower Partial Moment (FLPM) (2)",
        "Second Lower Partial Moment (SLPM) (2)",
        "Value at Risk (VaR) (2)",
        "Conditional Value at Risk (CVaR) (2)",
        "Entropic Value at Risk (EVaR) (2)",
        "Tail Gini of Losses (TG) (2)",
        "Relativistic Value at Risk (RLVaR ) (2)",
        "Worst Realization (2)",
        "Skewness",
        "Kurtosis",
        "",
        "Risk Measures based on Drawdowns (3)",
        "Ulcer Index (UCI)",
        "Average Drawdown (ADD)",
        "Drawdown at Risk (DaR)",
        "Conditional Drawdown at Risk (CDaR)",
        "Entropic Drawdown at Risk (EDaR)",
        "Relativistic Drawdown at Risk (RLDaR)",
        "Max Drawdown (MDD)",
        "(1) Annualized, multiplied by " + str(t_factor),
        "(2) Annualized, multiplied by √" + str(t_factor),
        "(3) Based on uncompounded cumulated returns",
    ]

    indicators = [
        "",
        (mu @ w_).to_numpy().item() * t_factor,
        np.power(np.prod(1 + X), days_per_year / days) - 1,
        MAR * t_factor,
        alpha,
        "",
        "",
        np.sqrt(w_.T @ cov @ w_).to_numpy().item() * t_factor**0.5,
        rk.MAD(X) * t_factor**0.5,
        rk.SemiDeviation(X) * t_factor**0.5,
        rk.LPM(X, MAR=MAR, p=1) * t_factor**0.5,
        rk.LPM(X, MAR=MAR, p=2) * t_factor**0.5,
        rk.VaR_Hist(X, alpha=alpha) * t_factor**0.5,
        rk.CVaR_Hist(X, alpha=alpha) * t_factor**0.5,
        rk.TG(X, alpha=alpha, a_sim=a_sim) * t_factor**0.5,
        rk.EVaR_Hist(X, alpha=alpha)[0] * t_factor**0.5,
        rk.RLVaR_Hist(X, alpha=alpha, kappa=kappa, solver=solver) * t_factor**0.5,
        rk.WR(X) * t_factor**0.5,
        st.skew(X, bias=False),
        st.kurtosis(X, bias=False),
        "",
        "",
        rk.UCI_Abs(X),
        rk.ADD_Abs(X),
        rk.DaR_Abs(X),
        rk.CDaR_Abs(X, alpha=alpha),
        rk.EDaR_Abs(X, alpha=alpha)[0],
        rk.RLDaR_Abs(X, alpha=alpha, kappa=kappa, solver=solver),
        rk.MDD_Abs(X),
        "",
        "",
        "",
    ]

    ratios = []
    for i in range(len(indicators)):
        if i < 6 or indicators[i] == "" or rowLabels[i] in ["Skewness", "Kurtosis"]:
            ratios.append("")
        else:
            if indicators[i] == 0:
                ratios.append("")
            else:
                ratio = (indicators[1] - MAR) / indicators[i]
                ratios.append(ratio)

    for i in range(len(indicators)):
        if indicators[i] != "":
            if rowLabels[i] in ["Skewness", "Kurtosis"]:
                indicators[i] = "{:.5f}".format(indicators[i])
            else:
                indicators[i] = "{:.4%}".format(indicators[i])
        if ratios[i] != "":
            ratios[i] = "{:.6f}".format(ratios[i])

    data = pd.DataFrame({"A": rowLabels, "B": indicators, "C": ratios}).to_numpy()

    ax.set_axis_off()
    ax.axis("tight")
    ax.axis("off")

    colLabels = ["", "Values", "(Return - MAR)/Risk"]
    colWidths = [0.45, 0.275, 0.275]
    rowHeight = 0.07

    table = ax.table(
        cellText=data,
        colLabels=colLabels,
        colWidths=colWidths,
        cellLoc="center",
        loc="upper left",
        bbox=[-0.03, 0, 1, 1],
    )

    table.auto_set_font_size(False)

    cellDict = table.get_celld()
    k = 1

    rowHeight = 1 / len(rowLabels)
    ncols = len(colLabels)
    nrows = len(rowLabels)

    for i in range(0, ncols):
        cellDict[(0, i)].set_text_props(weight="bold", color="white", size="x-large")
        cellDict[(0, i)].set_facecolor("darkblue")
        cellDict[(0, i)].set_edgecolor("white")
        cellDict[(0, i)].set_height(rowHeight)
        for j in range(1, nrows + 1):
            cellDict[(j, 0)].set_text_props(
                weight="bold", color="black", size="x-large", ha="left"
            )
            cellDict[(j, i)].set_text_props(color="black", size="x-large")
            cellDict[(j, 0)].set_edgecolor("white")
            cellDict[(j, i)].set_edgecolor("white")
            if k % 2 != 0:
                cellDict[(j, 0)].set_facecolor("whitesmoke")
                cellDict[(j, i)].set_facecolor("whitesmoke")
            if j in [6, 20]:
                cellDict[(j, 0)].set_facecolor("white")
                cellDict[(j, i)].set_facecolor("white")
            if j in [1, 7, 21]:
                cellDict[(j, 0)].set_text_props(color="white")
                cellDict[(j, 0)].set_facecolor("orange")
                cellDict[(j, i)].set_facecolor("orange")
                k = 1
            k += 1

            cellDict[(j, i)].set_height(rowHeight)

    for i in range(0, ncols):
        for j in range(nrows - 2, nrows + 1):
            cellDict[(j, i)].set_text_props(
                weight="normal", color="black", size="large"
            )
            cellDict[(j, i)].set_facecolor("white")

    try:
        fig.set_layout_engine(layout="constrained")
    except:
        pass

    return ax


def plot_clusters(
    returns,
    custom_cov=None,
    codependence="pearson",
    linkage="ward",
    opt_k_method="twodiff",
    k=None,
    max_k=10,
    bins_info="KN",
    alpha_tail=0.05,
    gs_threshold=0.5,
    leaf_order=True,
    show_clusters=True,
    dendrogram=True,
    cmap="RdYlBu",
    linecolor="fuchsia",
    title="",
    height=11,
    width=12,
    ax=None,
):
    r"""
    Create a clustermap plot based on the selected codependence measure.

    Parameters
    ----------
    returns : DataFrame of shape (n_samples, n_assets)
        Assets returns DataFrame, where n_samples is the number of
        observations and n_assets is the number of assets.
    custom_cov : DataFrame or None, optional
        Custom covariance matrix, used when codependence parameter has value
        'custom_cov'. The default is None.
    codependence : str, can be {'pearson', 'spearman', 'abs_pearson', 'abs_spearman', 'distance', 'mutual_info', 'tail' or 'custom_cov'}
        The codependence or similarity matrix used to build the distance
        metric and clusters. The default is 'pearson'. Possible values are:

        - 'pearson': pearson correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{0.5(1-\rho_{i,j})}`.
        - 'spearman': spearman correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{0.5(1-\rho_{i,j})}`.
        - 'kendall': kendall correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{0.5(1-\rho^{kendall}_{i,j})}`.
        - 'gerber1': Gerber statistic 1 correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{0.5(1-\rho^{gerber1}_{i,j})}`.
        - 'gerber2': Gerber statistic 2 correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{0.5(1-\rho^{gerber2}_{i,j})}`.
        - 'abs_pearson': absolute value pearson correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{(1-|\rho_{i,j}|)}`.
        - 'abs_spearman': absolute value spearman correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{(1-|\rho_{i,j}|)}`.
        - 'abs_kendall': absolute value kendall correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{(1-|\rho^{kendall}_{i,j}|)}`.
        - 'distance': distance correlation matrix. Distance formula :math:`D_{i,j} = \sqrt{(1-|\rho_{i,j}|)}`.
        - 'mutual_info': mutual information matrix. Distance used is variation information matrix.
        - 'tail': lower tail dependence index matrix. Dissimilarity formula :math:`D_{i,j} = -\log{\lambda_{i,j}}`.
        - 'custom_cov': use custom correlation matrix based on the custom_cov parameter. Distance formula: :math:`D_{i,j} = \sqrt{0.5(1-\rho^{pearson}_{i,j})}`.

    linkage : string, optional
        Linkage method of hierarchical clustering, see `linkage <https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html?highlight=linkage#scipy.cluster.hierarchy.linkage>`_ for more details.
        The default is 'ward'. Possible values are:

        - 'single'.
        - 'complete'.
        - 'average'.
        - 'weighted'.
        - 'centroid'.
        - 'median'.
        - 'ward'.
        - 'DBHT': Direct Bubble Hierarchical Tree.

    opt_k_method : str
        Method used to calculate the optimum number of clusters.
        The default is 'twodiff'. Possible values are:

        - 'twodiff': two difference gap statistic.
        - 'stdsil': standarized silhouette score.

    k : int, optional
        Number of clusters. This value is took instead of the optimal number
        of clusters calculated with the two difference gap statistic.
        The default is None.
    max_k : int, optional
        Max number of clusters used by the two difference gap statistic
        to find the optimal number of clusters. The default is 10.
    bins_info: int or str
        Number of bins used to calculate variation of information. The default
        value is 'KN'. Possible values are:

        - 'KN': Knuth's choice method. See more in `knuth_bin_width <https://docs.astropy.org/en/stable/api/astropy.stats.knuth_bin_width.html>`_.
        - 'FD': Freedman–Diaconis' choice method. See more in `freedman_bin_width <https://docs.astropy.org/en/stable/api/astropy.stats.freedman_bin_width.html>`_.
        - 'SC': Scotts' choice method. See more in `scott_bin_width <https://docs.astropy.org/en/stable/api/astropy.stats.scott_bin_width.html>`_.
        - 'HGR': Hacine-Gharbi and Ravier' choice method.
        - int: integer value chosen by user.

    alpha_tail : float, optional
        Significance level for lower tail dependence index. The default is 0.05.
    gs_threshold : float, optional
        Gerber statistic threshold. The default is 0.5.
    leaf_order : bool, optional
        Indicates if the cluster are ordered so that the distance between
        successive leaves is minimal. The default is True.
    show_clusters : bool, optional
        Indicates if clusters are plot. The default is True.
    dendrogram : bool, optional
        Indicates if the plot has or not a dendrogram. The default is True.
    cmap : str or cmap, optional
        Colormap used to plot the pcolormesh plot. The default is 'viridis'.
    linecolor : str, optional
        Color used to identify the clusters in the pcolormesh plot.
        The default is fuchsia'.
    title : str, optional
        Title of the chart. The default is "".
    height : float, optional
        Height of the image in inches. The default is 12.
    width : float, optional
        Width of the image in inches. The default is 12.
    ax : matplotlib axis, optional
        If provided, plot on this axis. The default is None.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    Returns
    -------
    ax : matplotlib axis
        Returns the Axes object with the plot for further tweaking.

    Example
    -------
    ::

        ax = rp.plot_clusters(returns=Y,
                              codependence='spearman',
                              linkage='ward',
                              k=None,
                              max_k=10,
                              leaf_order=True,
                              dendrogram=True,
                              ax=None)

    .. image:: images/Assets_Clusters.png


    """

    if not isinstance(returns, pd.DataFrame):
        raise ValueError("returns must be a DataFrame")

    if custom_cov is not None:
        if isinstance(custom_cov, pd.DataFrame):
            if custom_cov.index.tolist() != custom_cov.columns.tolist():
                raise ValueError(
                    "custom_cov must be a square DataFrame with samen labels in indexes and columns."
                )
            elif returns.columns.tolist() != custom_cov.index.tolist():
                raise ValueError("returns and custom_cov must have the same assets.")
        else:
            raise ValueError("custom_cov must be a square DataFrame.")

    width_ratios_1 = [0.25, 0.7, 0.05]
    height_ratios_1 = [0.27, 0.73]
    width_ratios_2 = [0.7, 0.05]

    if ax is None:
        fig = plt.gcf()
        ax = fig.gca()
        ax.axis("off")
        if dendrogram == True:
            gs = GridSpec(
                nrows=2,
                ncols=3,
                figure=fig,
                hspace=0.01,
                wspace=0.01,
                height_ratios=height_ratios_1,
                width_ratios=width_ratios_1,
            )
            axes = []
            for i in range(2):
                for j in range(3):
                    axes.append(fig.add_subplot(gs[i, j]))
        else:
            gs = GridSpec(nrows=1, ncols=2, figure=fig, width_ratios=width_ratios_2)
            axes = []
            axes.append(fig.add_subplot(gs[0, 0]))
            axes.append(fig.add_subplot(gs[0, 1]))
        fig.set_figwidth(width)
        fig.set_figheight(height)
    else:
        ax.axis("off")
        fig = ax.get_figure()
        if dendrogram == True:
            if isinstance(ax, plt.Axes):
                if hasattr(ax, "get_subplotspec"):
                    subplot_spec = ax.get_subplotspec()
                    gs0 = subplot_spec.get_gridspec()
                    gs = GridSpecFromSubplotSpec(
                        nrows=2,
                        ncols=3,
                        height_ratios=height_ratios_1,
                        width_ratios=width_ratios_1,
                        subplot_spec=gs0[0],
                    )
                else:
                    gs = GridSpec(
                        nrows=2,
                        ncols=3,
                        figure=fig,
                        height_ratios=height_ratios_1,
                        width_ratios=width_ratios_1,
                    )
                axes = []
                for i in range(2):
                    for j in range(3):
                        axes.append(fig.add_subplot(gs[i, j]))
            else:
                raise TypeError("ax must be a matplotlib axes object.")
        else:
            if isinstance(ax, plt.Axes):
                if hasattr(ax, "get_subplotspec"):
                    subplot_spec = ax.get_subplotspec()
                    gs0 = subplot_spec.get_gridspec()
                    gs = GridSpecFromSubplotSpec(
                        nrows=1,
                        ncols=2,
                        width_ratios=width_ratios_2,
                        subplot_spec=gs0[0],
                    )
                else:
                    gs = GridSpec(
                        nrows=1,
                        ncols=2,
                        figure=fig,
                        width_ratios=width_ratios_2,
                    )
                axes = []
                for i in range(1):
                    for j in range(2):
                        axes.append(fig.add_subplot(gs[i, j]))
            else:
                raise TypeError("ax must be a matplotlib axes object.")

    for i in range(len(axes) - 1):
        axes[i].grid(False)
        axes[i].axis("off")

    if dendrogram == True:
        (
            ax0,
            axcolor,
        ) = (
            axes[4],
            axes[5],
        )
    else:
        (
            ax0,
            axcolor,
        ) = (
            axes[0],
            axes[1],
        )

    labels = np.array(returns.columns.tolist())

    vmin, vmax = 0, 1
    if codependence in {
        "pearson",
        "spearman",
        "kendall",
        "gerber1",
        "gerber2",
    }:
        vmin, vmax = -1, 1

    # Calculating codependence matrix and distance metric
    codep, dist = af.codep_dist(
        returns=returns,
        custom_cov=custom_cov,
        codependence=codependence,
        bins_info=bins_info,
        alpha_tail=alpha_tail,
        gs_threshold=gs_threshold,
    )

    # Hierarchical clustering
    dist = dist.to_numpy()
    dist = pd.DataFrame(dist, columns=codep.columns, index=codep.index)
    dim = len(dist)
    if linkage == "DBHT":
        # different choices for D, S give different outputs!
        D = dist.to_numpy()  # dissimilarity matrix
        if codependence in {"pearson", "spearman", "custom_cov"}:
            S = (1 - dist**2).to_numpy()
        else:
            S = codep.to_numpy()  # similarity matrix
        (_, _, _, _, _, clustering) = db.DBHTs(
            D, S, leaf_order=leaf_order
        )  # DBHT clustering
    else:
        p_dist = squareform(dist, checks=False)
        clustering = hr.linkage(p_dist, method=linkage, optimal_ordering=leaf_order)

    # Ordering clusterings
    permutation = hr.leaves_list(clustering)
    permutation = permutation.tolist()
    ordered_codep = codep.to_numpy()[permutation, :][:, permutation]

    # optimal number of clusters
    if k is None:
        if opt_k_method == "twodiff":
            k, clustering_inds = af.two_diff_gap_stat(dist, clustering, max_k)
        elif opt_k_method == "stdsil":
            k, clustering_inds = af.std_silhouette_score(dist, clustering, max_k)
        else:
            raise ValueError("The only opt_k_method available are twodiff and stdsil")
    else:
        clustering_inds = hr.fcluster(clustering, k, criterion="maxclust")

    clusters = {i: [] for i in np.unique(clustering_inds)}
    for i, v in enumerate(clustering_inds):
        clusters[v].append(i)

    # ax[4] = fig.add_axes([0.2, 0.15, 0.55, 0.55])

    im = ax0.pcolormesh(ordered_codep, cmap=cmap, vmin=vmin, vmax=vmax)
    ax0.axis("on")
    ax0.set_xticks(np.arange(codep.shape[0]) + 0.5, minor=False)
    ax0.set_yticks(np.arange(codep.shape[0]) + 0.5, minor=False)
    ax0.set_xticklabels(labels[permutation], rotation=90, ha="center")
    ax0.set_yticklabels(labels[permutation], va="center")
    ax0.yaxis.set_label_position("right")
    ax0.yaxis.tick_right()
    ax0.set_ylim(ax0.get_ylim()[::-1])

    flag = False
    if show_clusters is True:
        if linecolor is None:
            linecolor = "fuchsia"
            flag = True
        elif linecolor is not None:
            flag = True

    if flag:
        N = len(permutation)
        for cluster_id, cluster in clusters.items():
            amin = permutation.index(cluster[0])
            xmin, xmax = amin, amin + len(cluster)
            ymin, ymax = amin, amin + len(cluster)

            for i in cluster:
                a = permutation.index(i)
                if a < amin:
                    xmin, xmax = a, a + len(cluster)
                    ymin, ymax = a, a + len(cluster)
                    amin = a

            ax0.axvline(
                x=xmin,
                ymin=(N - ymin) / dim,
                ymax=(N - ymax) / dim,
                linewidth=4,
                color=linecolor,
            )
            ax0.axvline(
                x=xmax,
                ymin=(N - ymin) / dim,
                ymax=(N - ymax) / dim,
                linewidth=4,
                color=linecolor,
            )
            ax0.axhline(
                y=ymin, xmin=xmin / dim, xmax=xmax / dim, linewidth=4, color=linecolor
            )
            ax0.axhline(
                y=ymax, xmin=xmin / dim, xmax=xmax / dim, linewidth=4, color=linecolor
            )

    # axcolor = fig.add_axes([.87, 0.15, 0.02, 0.55])
    ax0.get_figure().colorbar(im, cax=axcolor)

    if dendrogram == True:
        # ax1 = fig.add_axes([0.2, 0.71, 0.55, 0.2])
        ax1 = axes[1]
        if show_clusters is False:
            color_threshold = 0
        elif show_clusters is True:
            L, M = hr.leaders(clustering, np.array(clustering_inds, dtype=np.int32))
            root, nodes = hr.to_tree(clustering, rd=True)
            nodes = np.array([i.dist for i in nodes])
            nodes.sort()
            leaders_threshold = nodes[np.max(L) + 1]
            color_threshold = np.max(leaders_threshold)
            colors = af.color_list(k)
            hr.set_link_color_palette(colors)

        hr.dendrogram(
            clustering,
            color_threshold=color_threshold,
            above_threshold_color="grey",
            ax=ax1,
        )
        hr.set_link_color_palette(None)
        ax1.xaxis.set_major_locator(mticker.FixedLocator(np.arange(codep.shape[0])))
        ax1.set_xticklabels(labels[permutation], rotation=90, ha="center")

        if show_clusters is True:
            i = 0
            for coll in ax1.collections[
                :-1
            ]:  # the last collection is the ungrouped level
                xmin, xmax = np.inf, -np.inf
                ymax = -np.inf
                for p in coll.get_paths():
                    (x0, _), (x1, y1) = p.get_extents().get_points()
                    xmin = min(xmin, x0)
                    xmax = max(xmax, x1)
                    ymax = max(ymax, y1)
                rec = plt.Rectangle(
                    (xmin - 4, 0),
                    xmax - xmin + 8,
                    ymax * 1.005,
                    facecolor=colors[i],  # coll.get_color()[0],
                    alpha=0.2,
                    edgecolor="none",
                )
                ax1.add_patch(rec)
                i += 1

        ax1.set_xticks([])
        ax1.set_yticks([])

        for i in {"right", "left", "top", "bottom"}:
            side = ax1.spines[i]
            side.set_visible(False)

        # ax2 = fig.add_axes([0.0, 0.15, 0.2, 0.55])
        ax2 = axes[3]
        if show_clusters is True:
            hr.set_link_color_palette(colors)

        hr.dendrogram(
            clustering,
            color_threshold=color_threshold,
            above_threshold_color="grey",
            orientation="left",
            ax=ax2,
        )
        hr.set_link_color_palette(None)

        ax2.xaxis.set_major_locator(mticker.FixedLocator(np.arange(codep.shape[0])))
        ax2.set_xticklabels(labels[permutation], rotation=90, ha="center")
        ax2.set_ylim(ax2.get_ylim()[::-1])

        if show_clusters is True:
            i = 0
            # the last collection is the ungrouped level
            for coll in ax2.collections[:-1]:
                ymin, ymax = np.inf, -np.inf
                xmax = -np.inf
                for p in coll.get_paths():
                    (_, y0), (x1, y1) = p.get_extents().get_points()
                    ymin = min(ymin, y0)
                    ymax = max(ymax, y1)
                    xmax = max(xmax, x1)
                rec = plt.Rectangle(
                    (0, ymin - 4),
                    xmax * 1.05,
                    ymax - ymin + 8,
                    facecolor=colors[i],  # coll.get_color()[0],
                    alpha=0.2,
                    edgecolor="none",
                )
                ax2.add_patch(rec)
                i += 1

        ax2.set_xticks([])
        ax2.set_yticks([])
        ax2.set_yticklabels([])
        for i in {"right", "left", "top", "bottom"}:
            side = ax2.spines[i]
            side.set_visible(False)

    if title == "":
        title = (
            "Assets Clustermap ("
            + codependence.capitalize()
            + " & "
            + linkage
            + " linkage)"
        )

    if dendrogram == True:
        ax1.set_title(title)
    elif dendrogram == False:
        ax0.set_title(title)

    try:
        fig.set_layout_engine(layout="constrained")
    except:
        pass

    return ax


def plot_dendrogram(
    returns,
    custom_cov=None,
    codependence="pearson",
    linkage="ward",
    opt_k_method="twodiff",
    k=None,
    max_k=10,
    bins_info="KN",
    alpha_tail=0.05,
    gs_threshold=0.5,
    leaf_order=True,
    show_clusters=True,
    title="",
    height=5,
    width=12,
    ax=None,
):
    r"""
    Create a dendrogram based on the selected codependence measure.

    Parameters
    ----------
    returns : DataFrame of shape (n_samples, n_assets)
        Assets returns DataFrame, where n_samples is the number of
        observations and n_assets is the number of assets.
    custom_cov : DataFrame or None, optional
        Custom covariance matrix, used when codependence parameter has value
        'custom_cov'. The default is None.
    codependence : str, can be {'pearson', 'spearman', 'abs_pearson', 'abs_spearman', 'distance', 'mutual_info', 'tail' or 'custom_cov'}
        The codependence or similarity matrix used to build the distance
        metric and clusters. The default is 'pearson'. Possible values are:

        - 'pearson': pearson correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{0.5(1-\rho_{i,j})}`.
        - 'spearman': spearman correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{0.5(1-\rho_{i,j})}`.
        - 'kendall': kendall correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{0.5(1-\rho^{kendall}_{i,j})}`.
        - 'gerber1': Gerber statistic 1 correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{0.5(1-\rho^{gerber1}_{i,j})}`.
        - 'gerber2': Gerber statistic 2 correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{0.5(1-\rho^{gerber2}_{i,j})}`.
        - 'abs_pearson': absolute value pearson correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{(1-|\rho_{i,j}|)}`.
        - 'abs_spearman': absolute value spearman correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{(1-|\rho_{i,j}|)}`.
        - 'abs_kendall': absolute value kendall correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{(1-|\rho^{kendall}_{i,j}|)}`.
        - 'distance': distance correlation matrix. Distance formula :math:`D_{i,j} = \sqrt{(1-|\rho_{i,j}|)}`.
        - 'mutual_info': mutual information matrix. Distance used is variation information matrix.
        - 'tail': lower tail dependence index matrix. Dissimilarity formula :math:`D_{i,j} = -\log{\lambda_{i,j}}`.
        - 'custom_cov': use custom correlation matrix based on the custom_cov parameter. Distance formula: :math:`D_{i,j} = \sqrt{0.5(1-\rho^{pearson}_{i,j})}`.

    linkage : string, optional
        Linkage method of hierarchical clustering, see `linkage <https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html?highlight=linkage#scipy.cluster.hierarchy.linkage>`_ for more details.
        The default is 'ward'. Possible values are:

        - 'single'.
        - 'complete'.
        - 'average'.
        - 'weighted'.
        - 'centroid'.
        - 'median'.
        - 'ward'.
        - 'DBHT': Direct Bubble Hierarchical Tree.

    opt_k_method : str
        Method used to calculate the optimum number of clusters.
        The default is 'twodiff'. Possible values are:

        - 'twodiff': two difference gap statistic.
        - 'stdsil': standarized silhouette score.

    k : int, optional
        Number of clusters. This value is took instead of the optimal number
        of clusters calculated with the two difference gap statistic.
        The default is None.
    max_k : int, optional
        Max number of clusters used by the two difference gap statistic
        to find the optimal number of clusters. The default is 10.
    bins_info: int or str
        Number of bins used to calculate variation of information. The default
        value is 'KN'. Possible values are:

        - 'KN': Knuth's choice method. See more in `knuth_bin_width <https://docs.astropy.org/en/stable/api/astropy.stats.knuth_bin_width.html>`_.
        - 'FD': Freedman–Diaconis' choice method. See more in `freedman_bin_width <https://docs.astropy.org/en/stable/api/astropy.stats.freedman_bin_width.html>`_.
        - 'SC': Scotts' choice method. See more in `scott_bin_width <https://docs.astropy.org/en/stable/api/astropy.stats.scott_bin_width.html>`_.
        - 'HGR': Hacine-Gharbi and Ravier' choice method.
        - int: integer value chosen by user.

    alpha_tail : float, optional
        Significance level for lower tail dependence index. The default is 0.05.
    gs_threshold : float, optional
        Gerber statistic threshold. The default is 0.5.
    leaf_order : bool, optional
        Indicates if the cluster are ordered so that the distance between
        successive leaves is minimal. The default is True.
    show_clusters : bool, optional
        Indicates if clusters are plot. The default is True.
    title : str, optional
        Title of the chart. The default is "".
    height : float, optional
        Height of the image in inches. The default is 5.
    width : float, optional
        Width of the image in inches. The default is 12.
    ax : matplotlib axis, optional
        If provided, plot on this axis. The default is None.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    Returns
    -------
    ax : matplotlib axis
        Returns the Axes object with the plot for further tweaking.

    Example
    -------
    ::

        ax = rp.plot_dendrogram(returns=Y,
                                codependence='spearman',
                                linkage='ward',
                                k=None,
                                max_k=10,
                                leaf_order=True,
                                ax=None)

    .. image:: images/Assets_Dendrogram.png


    """
    if not isinstance(returns, pd.DataFrame):
        raise ValueError("returns must be a DataFrame")

    if custom_cov is not None:
        if isinstance(custom_cov, pd.DataFrame):
            if custom_cov.index.tolist() != custom_cov.columns.tolist():
                raise ValueError(
                    "custom_cov must be a square DataFrame with samen labels in indexes and columns."
                )
            elif returns.columns.tolist() != custom_cov.index.tolist():
                raise ValueError("returns and custom_cov must have the same assets.")
        else:
            raise ValueError("custom_cov must be a square DataFrame.")

    if ax is None:
        fig = plt.gcf()
        ax = fig.gca()
        fig.set_figwidth(width)
        fig.set_figheight(height)
    else:
        fig = ax.get_figure()

    labels = np.array(returns.columns.tolist())

    # Calculating codependence matrix and distance metric
    codep, dist = af.codep_dist(
        returns=returns,
        custom_cov=custom_cov,
        codependence=codependence,
        bins_info=bins_info,
        alpha_tail=alpha_tail,
        gs_threshold=gs_threshold,
    )

    # Hierarchical clustering
    dist = dist.to_numpy()
    dist = pd.DataFrame(dist, columns=codep.columns, index=codep.index)
    if linkage == "DBHT":
        # different choices for D, S give different outputs!
        D = dist.to_numpy()  # dissimilarity matrix
        if codependence in {"pearson", "spearman", "custom_cov"}:
            S = (1 - dist**2).to_numpy()
        else:
            S = codep.copy().to_numpy()  # similarity matrix
        (_, _, _, _, _, clustering) = db.DBHTs(
            D, S, leaf_order=leaf_order
        )  # DBHT clustering
    else:
        p_dist = squareform(dist, checks=False)
        clustering = hr.linkage(p_dist, method=linkage, optimal_ordering=leaf_order)

    # Ordering clusterings
    permutation = hr.leaves_list(clustering)
    permutation = permutation.tolist()

    if show_clusters is False:
        color_threshold = 0
    elif show_clusters is True:
        # optimal number of clusters
        if k is None:
            if opt_k_method == "twodiff":
                k, clustering_inds = af.two_diff_gap_stat(dist, clustering, max_k)
            elif opt_k_method == "stdsil":
                k, clustering_inds = af.std_silhouette_score(dist, clustering, max_k)
            else:
                raise ValueError(
                    "The only opt_k_method available are twodiff and stdsil"
                )
        else:
            clustering_inds = hr.fcluster(clustering, k, criterion="maxclust")

        L, M = hr.leaders(clustering, np.array(clustering_inds, dtype=np.int32))
        root, nodes = hr.to_tree(clustering, rd=True)
        nodes = np.array([i.dist for i in nodes])
        nodes.sort()
        leaders_threshold = nodes[np.max(L) + 1]
        color_threshold = np.max(leaders_threshold)
        colors = af.color_list(k)  # color list
        hr.set_link_color_palette(colors)

    hr.dendrogram(
        clustering, color_threshold=color_threshold, above_threshold_color="grey", ax=ax
    )
    hr.set_link_color_palette(None)

    ax.set_xticklabels(labels[permutation], rotation=90, ha="center")

    if show_clusters is True:
        i = 0
        for coll in ax.collections[:-1]:  # the last collection is the ungrouped level
            xmin, xmax = np.inf, -np.inf
            ymax = -np.inf
            for p in coll.get_paths():
                (x0, _), (x1, y1) = p.get_extents().get_points()
                xmin = min(xmin, x0)
                xmax = max(xmax, x1)
                ymax = max(ymax, y1)
            rec = plt.Rectangle(
                (xmin - 4, 0),
                xmax - xmin + 8,
                ymax * 1.005,
                facecolor=colors[i],  # coll.get_color()[0],
                alpha=0.2,
                edgecolor="none",
            )
            ax.add_patch(rec)
            i += 1

    ax.set_yticks([])
    ax.set_yticklabels([])
    for i in {"right", "left", "top", "bottom"}:
        side = ax.spines[i]
        side.set_visible(False)

    if title == "":
        title = (
            "Assets Dendrogram ("
            + codependence.capitalize()
            + " & "
            + linkage
            + " linkage)"
        )

    ax.set_title(title)

    try:
        fig.set_layout_engine(layout="constrained")
    except:
        pass

    return ax


def plot_network(
    returns,
    custom_cov=None,
    codependence="pearson",
    linkage="ward",
    opt_k_method="twodiff",
    k=None,
    max_k=10,
    bins_info="KN",
    alpha_tail=0.05,
    gs_threshold=0.5,
    leaf_order=True,
    kind="spring",
    seed=0,
    node_labels=True,
    node_size=1400,
    node_alpha=0.7,
    font_size=10,
    title="",
    height=8,
    width=10,
    ax=None,
):
    r"""
    Create a network plot. The Planar Maximally Filtered Graph (PMFG) for DBHT
    linkage and Minimum Spanning Tree (MST) for other linkage methods.

    Parameters
    ----------
    returns : DataFrame of shape (n_samples, n_assets)
        Assets returns DataFrame, where n_samples is the number of
        observations and n_assets is the number of assets.
    custom_cov : DataFrame or None, optional
        Custom covariance matrix, used when codependence parameter has value
        'custom_cov'. The default is None.
    codependence : str, can be {'pearson', 'spearman', 'abs_pearson', 'abs_spearman', 'distance', 'mutual_info', 'tail' or 'custom_cov'}
        The codependence or similarity matrix used to build the distance
        metric and clusters. The default is 'pearson'. Possible values are:

        - 'pearson': pearson correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{0.5(1-\rho^{pearson}_{i,j})}`.
        - 'spearman': spearman correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{0.5(1-\rho^{spearman}_{i,j})}`.
        - 'kendall': kendall correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{0.5(1-\rho^{kendall}_{i,j})}`.
        - 'gerber1': Gerber statistic 1 correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{0.5(1-\rho^{gerber1}_{i,j})}`.
        - 'gerber2': Gerber statistic 2 correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{0.5(1-\rho^{gerber2}_{i,j})}`.
        - 'abs_pearson': absolute value pearson correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{(1-|\rho_{i,j}|)}`.
        - 'abs_spearman': absolute value spearman correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{(1-|\rho_{i,j}|)}`.
        - 'abs_kendall': absolute value kendall correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{(1-|\rho^{kendall}_{i,j}|)}`.
        - 'distance': distance correlation matrix. Distance formula :math:`D_{i,j} = \sqrt{(1-\rho^{distance}_{i,j})}`.
        - 'mutual_info': mutual information matrix. Distance used is variation information matrix.
        - 'tail': lower tail dependence index matrix. Dissimilarity formula :math:`D_{i,j} = -\log{\lambda_{i,j}}`.
        - 'custom_cov': use custom correlation matrix based on the custom_cov parameter. Distance formula: :math:`D_{i,j} = \sqrt{0.5(1-\rho^{pearson}_{i,j})}`.

    linkage : string, optional
        Linkage method of hierarchical clustering, see `linkage <https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html?highlight=linkage#scipy.cluster.hierarchy.linkage>`_ for more details.
        The default is 'ward'. Possible values are:

        - 'single'.
        - 'complete'.
        - 'average'.
        - 'weighted'.
        - 'centroid'.
        - 'median'.
        - 'ward'.
        - 'DBHT': Direct Bubble Hierarchical Tree.

    opt_k_method : str
        Method used to calculate the optimum number of clusters.
        The default is 'twodiff'. Possible values are:

        - 'twodiff': two difference gap statistic.
        - 'stdsil': standarized silhouette score.

    k : int, optional
        Number of clusters. This value is took instead of the optimal number
        of clusters calculated with the two difference gap statistic.
        The default is None.
    max_k : int, optional
        Max number of clusters used by the two difference gap statistic
        to find the optimal number of clusters. The default is 10.
    bins_info: int or str
        Number of bins used to calculate variation of information. The default
        value is 'KN'. Possible values are:

        - 'KN': Knuth's choice method. See more in `knuth_bin_width <https://docs.astropy.org/en/stable/api/astropy.stats.knuth_bin_width.html>`_.
        - 'FD': Freedman–Diaconis' choice method. See more in `freedman_bin_width <https://docs.astropy.org/en/stable/api/astropy.stats.freedman_bin_width.html>`_.
        - 'SC': Scotts' choice method. See more in `scott_bin_width <https://docs.astropy.org/en/stable/api/astropy.stats.scott_bin_width.html>`_.
        - 'HGR': Hacine-Gharbi and Ravier' choice method.
        - int: integer value chosen by user.

    alpha_tail : float, optional
        Significance level for lower tail dependence index. The default is 0.05.
    gs_threshold : float, optional
        Gerber statistic threshold. The default is 0.5.
    leaf_order : bool, optional
        Indicates if the cluster are ordered so that the distance between
        successive leaves is minimal. The default is True.
    kind : str, optional
        Kind of networkx layout. The default value is 'spring'. Possible values
        are:

        - 'spring': networkx spring_layout.
        - 'planar'. networkx planar_layout.
        - 'circular'. networkx circular_layout.
        - 'kamada'. networkx kamada_kawai_layout.

    seed : int, optional
        Seed for networkx spring layout. The default value is 0.
    node_labels : bool, optional
        Specify if node lables are visible. The default value is True.
    node_size : float, optional
        Size of the nodes. The default value is 1600.
    node_alpha : float, optional
        Alpha parameter or transparency of nodes. The default value is 0.7.
    font_size : float, optional
        Font size of node labels. The default value is 12.
    title : str, optional
        Title of the chart. The default is "".
    height : float, optional
        Height of the image in inches. The default is 8.
    width : float, optional
        Width of the image in inches. The default is 10.
    ax : matplotlib axis, optional
        If provided, plot on this axis. The default is None.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    Returns
    -------
    ax : matplotlib axis
        Returns the Axes object with the plot for further tweaking.

    Example
    -------
    ::

        ax = rp.plot_network(returns=Y,
                             codependence="pearson",
                             linkage="ward",
                             k=None,
                             max_k=10,
                             alpha_tail=0.05,
                             leaf_order=True,
                             kind='kamada',
                             ax=None)

    .. image:: images/Assets_Network.png


    """
    if not isinstance(returns, pd.DataFrame):
        raise ValueError("returns must be a DataFrame")

    if custom_cov is not None:
        if isinstance(custom_cov, pd.DataFrame):
            if custom_cov.index.tolist() != custom_cov.columns.tolist():
                raise ValueError(
                    "custom_cov must be a square DataFrame with samen labels in indexes and columns."
                )
            elif returns.columns.tolist() != custom_cov.index.tolist():
                raise ValueError("returns and custom_cov must have the same assets.")
        else:
            raise ValueError("custom_cov must be a square DataFrame.")

    if ax is None:
        fig = plt.gcf()
        ax = fig.gca()
        fig.set_figwidth(width)
        fig.set_figheight(height)
    else:
        fig = ax.get_figure()

    labels = np.array(returns.columns.tolist())

    # Calculating codependence matrix and distance metric
    codep, dist = af.codep_dist(
        returns=returns,
        custom_cov=custom_cov,
        codependence=codependence,
        bins_info=bins_info,
        alpha_tail=alpha_tail,
        gs_threshold=gs_threshold,
    )

    # Hierarchical clustering
    dist = dist.to_numpy()
    dist = pd.DataFrame(dist, columns=codep.columns, index=codep.index)
    if linkage == "DBHT":
        # different choices for D, S give different outputs!
        D = dist.to_numpy()  # dissimilarity matrix
        if codependence in {"pearson", "spearman", "custom_cov"}:
            S = (1 - dist**2).to_numpy()
        else:
            S = codep.copy().to_numpy()  # similarity matrix
        (_, Rpm, _, _, _, clustering) = db.DBHTs(
            D, S, leaf_order=leaf_order
        )  # DBHT clustering
        MAdj = pd.DataFrame(Rpm, index=labels, columns=labels)
        G = nx.from_pandas_adjacency(MAdj)
    else:
        p_dist = squareform(dist, checks=False)
        clustering = hr.linkage(p_dist, method=linkage, optimal_ordering=leaf_order)
        T = nx.from_pandas_adjacency(dist)  # create a graph G from a numpy matrix
        G = nx.minimum_spanning_tree(T)

    # optimal number of clusters
    if k is None:
        if opt_k_method == "twodiff":
            k, clustering_inds = af.two_diff_gap_stat(dist, clustering, max_k)
        elif opt_k_method == "stdsil":
            k, clustering_inds = af.std_silhouette_score(dist, clustering, max_k)
        else:
            raise ValueError("The only opt_k_method available are twodiff and stdsil")
    else:
        clustering_inds = hr.fcluster(clustering, k, criterion="maxclust")

    clusters = {i: [] for i in np.unique(clustering_inds)}
    for i, v in enumerate(clustering_inds):
        clusters[v].append(labels[i])

    # Layout options
    node_options = {
        "node_size": node_size,
        "alpha": node_alpha,
    }
    font_options = {
        "font_size": font_size,
        "font_color": "k",
    }

    label_options = {"ec": "k", "fc": "white", "alpha": 0.7}

    if kind == "spring":
        pos = nx.spring_layout(G, seed=int(seed))
    elif kind == "planar":
        pos = nx.planar_layout(G)
    elif kind == "circular":
        pos = nx.circular_layout(G)
    elif kind == "kamada":
        pos = nx.kamada_kawai_layout(G)

    # Plotting
    nx.draw_networkx_edges(G, pos=pos, ax=ax, edge_color="grey")

    if node_labels == True:
        nx.draw_networkx_labels(G, pos=pos, ax=ax, bbox=label_options, **font_options)

    colors = af.color_list(k)

    for i, color in zip(clusters.keys(), colors):
        nx.draw_networkx_nodes(
            G, pos=pos, nodelist=clusters[i], node_color=color, ax=ax, **node_options
        )

    ax.set_yticks([])
    ax.set_yticklabels([])
    for i in {"right", "left", "top", "bottom"}:
        side = ax.spines[i]
        side.set_visible(False)

    if title == "":
        if linkage == "DBHT":
            title = (
                "Planar Maximally Filtered Graph ("
                + codependence.capitalize()
                + ", "
                + linkage
                + " linkage & "
                + kind
                + " layout)"
            )
        else:
            title = (
                "Minimun Spanning Tree ("
                + codependence.capitalize()
                + ", "
                + linkage
                + " linkage & "
                + kind
                + " layout)"
            )

    ax.set_title(title)

    try:
        fig.set_layout_engine(layout="constrained")
    except:
        pass

    return ax


def plot_network_allocation(
    returns,
    w,
    custom_cov=None,
    codependence="pearson",
    linkage="ward",
    bins_info="KN",
    alpha_tail=0.05,
    gs_threshold=0.5,
    leaf_order=True,
    kind="spring",
    seed=0,
    node_labels=True,
    max_node_size=2000,
    color_lng="tab:blue",
    color_sht="tab:red",
    label_v=0.08,
    label_h=0,
    font_size=10,
    title="",
    height=8,
    width=10,
    ax=None,
):
    r"""
    Create a network plot with node size of the nodes and color represents the
    amount invested and direction (long-short) respectively. The Planar Maximally
    Filtered Graph (PMFG) for DBHT linkage and Minimum Spanning Tree (MST)
    for other linkage methods.

    Parameters
    ----------
    returns : DataFrame of shape (n_samples, n_assets)
        Assets returns DataFrame, where n_samples is the number of
        observations and n_assets is the number of assets.
    w : DataFrame or Series of shape (n_assets, 1)
        Portfolio weights, where n_assets is the number of assets.
    custom_cov : DataFrame or None, optional
        Custom covariance matrix, used when codependence parameter has value
        'custom_cov'. The default is None.
    codependence : str, can be {'pearson', 'spearman', 'abs_pearson', 'abs_spearman', 'distance', 'mutual_info', 'tail' or 'custom_cov'}
        The codependence or similarity matrix used to build the distance
        metric and clusters. The default is 'pearson'. Possible values are:

        - 'pearson': pearson correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{0.5(1-\rho^{pearson}_{i,j})}`.
        - 'spearman': spearman correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{0.5(1-\rho^{spearman}_{i,j})}`.
        - 'kendall': kendall correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{0.5(1-\rho^{kendall}_{i,j})}`.
        - 'gerber1': Gerber statistic 1 correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{0.5(1-\rho^{gerber1}_{i,j})}`.
        - 'gerber2': Gerber statistic 2 correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{0.5(1-\rho^{gerber2}_{i,j})}`.
        - 'abs_pearson': absolute value pearson correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{(1-|\rho_{i,j}|)}`.
        - 'abs_spearman': absolute value spearman correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{(1-|\rho_{i,j}|)}`.
        - 'abs_kendall': absolute value kendall correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{(1-|\rho^{kendall}_{i,j}|)}`.
        - 'distance': distance correlation matrix. Distance formula :math:`D_{i,j} = \sqrt{(1-\rho^{distance}_{i,j})}`.
        - 'mutual_info': mutual information matrix. Distance used is variation information matrix.
        - 'tail': lower tail dependence index matrix. Dissimilarity formula :math:`D_{i,j} = -\log{\lambda_{i,j}}`.
        - 'custom_cov': use custom correlation matrix based on the custom_cov parameter. Distance formula: :math:`D_{i,j} = \sqrt{0.5(1-\rho^{pearson}_{i,j})}`.

    linkage : string, optional
        Linkage method of hierarchical clustering, see `linkage <https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html?highlight=linkage#scipy.cluster.hierarchy.linkage>`_ for more details.
        The default is 'ward'. Possible values are:

        - 'single'.
        - 'complete'.
        - 'average'.
        - 'weighted'.
        - 'centroid'.
        - 'median'.
        - 'ward'.
        - 'DBHT': Direct Bubble Hierarchical Tree.

    bins_info: int or str
        Number of bins used to calculate variation of information. The default
        value is 'KN'. Possible values are:

        - 'KN': Knuth's choice method. See more in `knuth_bin_width <https://docs.astropy.org/en/stable/api/astropy.stats.knuth_bin_width.html>`_.
        - 'FD': Freedman–Diaconis' choice method. See more in `freedman_bin_width <https://docs.astropy.org/en/stable/api/astropy.stats.freedman_bin_width.html>`_.
        - 'SC': Scotts' choice method. See more in `scott_bin_width <https://docs.astropy.org/en/stable/api/astropy.stats.scott_bin_width.html>`_.
        - 'HGR': Hacine-Gharbi and Ravier' choice method.
        - int: integer value chosen by user.

    alpha_tail : float, optional
        Significance level for lower tail dependence index. The default is 0.05.
    gs_threshold : float, optional
        Gerber statistic threshold. The default is 0.5.
    leaf_order : bool, optional
        Indicates if the cluster are ordered so that the distance between
        successive leaves is minimal. The default is True.
    kind : str, optional
        Kind of networkx layout. The default value is 'spring'. Possible values
        are:

        - 'spring': networkx spring_layout.
        - 'planar'. networkx planar_layout.
        - 'circular'. networkx circular_layout.
        - 'kamada'. networkx kamada_kawai_layout.

    seed : int, optional
        Seed for networkx spring layout. The default value is 0.
    node_labels : bool, optional
        Specify if node lables are visible. The default value is True.
    max_node_size : float, optional
        Size of the node with maximum weight in absolute value. The default
        value is 2000.
    color_lng : str, optional
        Color of assets with long positions. The default value is 'tab:blue'.
    color_sht : str, optional
        Color of assets with short positions. The default value is 'tab:red'.
    label_v : float, optional
        Vertical distance the label is offset from the center. The default value is 0.08.
    label_h : float, optional
        Horizontal distance the label is offset from the center. The default value is 0.
    font_size : float, optional
        Font size of node labels. The default value is 12.
    title : str, optional
        Title of the chart. The default is "".
    height : float, optional
        Height of the image in inches. The default is 8.
    width : float, optional
        Width of the image in inches. The default is 10.
    ax : matplotlib axis, optional
        If provided, plot on this axis. The default is None.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    Returns
    -------
    ax : matplotlib axis
        Returns the Axes object with the plot for further tweaking.

    Example
    -------
    ::

        ax = rp.plot_network_allocation(returns=Y,
                                        w=w1,
                                        codependence="pearson",
                                        linkage="ward",
                                        alpha_tail=0.05,
                                        leaf_order=True,
                                        kind='kamada',
                                        ax=None)

    .. image:: images/Assets_Network_Allocation.png


    """
    if not isinstance(returns, pd.DataFrame):
        raise ValueError("returns must be a DataFrame")

    if not isinstance(w, pd.DataFrame):
        if isinstance(w, pd.Series):
            w_ = w.to_frame()
        else:
            raise ValueError("w must be a one column DataFrame or Series")
    else:
        if w.shape[0] == 1:
            w_ = w.T.copy()
        elif w.shape[1] == 1:
            w_ = w.copy()
        else:
            raise ValueError("w must be a one column DataFrame or Series")

    if returns.columns.tolist() != w_.index.tolist():
        raise ValueError("returns and w must have same assets.")

    if custom_cov is not None:
        if isinstance(custom_cov, pd.DataFrame):
            if custom_cov.index.tolist() != custom_cov.columns.tolist():
                raise ValueError(
                    "custom_cov must be a square DataFrame with samen labels in indexes and columns."
                )
            elif returns.columns.tolist() != custom_cov.index.tolist():
                raise ValueError("returns and custom_cov must have the same assets.")
        else:
            raise ValueError("custom_cov must be a square DataFrame.")

    if ax is None:
        fig = plt.gcf()
        ax = fig.gca()
        fig.set_figwidth(width)
        fig.set_figheight(height)
    else:
        fig = ax.get_figure()

    labels = np.array(returns.columns.tolist())

    # Calculating codependence matrix and distance metric
    codep, dist = af.codep_dist(
        returns=returns,
        custom_cov=custom_cov,
        codependence=codependence,
        bins_info=bins_info,
        alpha_tail=alpha_tail,
        gs_threshold=gs_threshold,
    )

    # Hierarchical clustering
    dist = dist.to_numpy()
    dist = pd.DataFrame(dist, columns=codep.columns, index=codep.index)
    if linkage == "DBHT":
        # different choices for D, S give different outputs!
        D = dist.to_numpy()  # dissimilarity matrix
        if codependence in {"pearson", "spearman", "custom_cov"}:
            S = (1 - dist**2).to_numpy()
        else:
            S = codep.copy().to_numpy()  # similarity matrix
        (_, Rpm, _, _, _, clustering) = db.DBHTs(
            D, S, leaf_order=leaf_order
        )  # DBHT clustering
        MAdj = pd.DataFrame(Rpm, index=labels, columns=labels)
        G = nx.from_pandas_adjacency(MAdj)
    else:
        T = nx.from_pandas_adjacency(dist)  # create a graph G from a numpy matrix
        G = nx.minimum_spanning_tree(T)

    label_options = {"ec": "k", "fc": "white", "alpha": 0.7}

    if kind == "spring":
        pos = nx.spring_layout(G, seed=int(seed))
    elif kind == "planar":
        pos = nx.planar_layout(G)
    elif kind == "circular":
        pos = nx.circular_layout(G)
    elif kind == "kamada":
        pos = nx.kamada_kawai_layout(G)

    # Plotting
    nx.draw_networkx_edges(G, pos=pos, ax=ax, edge_color="grey")

    clusters_1 = w_.loc[w_.iloc[:, 0] > 1e-6]
    clusters_2 = w_.loc[(w_.iloc[:, 0] <= 1e-6) & (w_.iloc[:, 0] >= -1e-6)]
    clusters_3 = w_.loc[w_.iloc[:, 0] < -1e-6]

    label_options = {"ec": "k", "fc": "white", "alpha": 0.7}

    font_options = {
        "font_size": font_size,
        "font_color": "k",
    }

    node_w = np.abs(w_) * max_node_size / np.abs(w_).sum().item()
    if node_labels == True:
        labels = {}
        labels_pos = {}
        for node in G.nodes():
            if node in list(clusters_1.index):
                labels[node] = node
                labels_pos[node] = pos[node] + np.array([label_h, label_v])

    nx.draw_networkx_nodes(
        G,
        pos=pos,
        nodelist=list(clusters_1.index),
        node_color=color_lng,
        alpha=0.8,
        ax=ax,
        node_size=node_w.loc[clusters_1.index],
    )
    nx.draw_networkx_nodes(
        G,
        pos=pos,
        nodelist=list(clusters_2.index),
        node_color="lightgrey",
        alpha=0.5,
        ax=ax,
        node_size=100,
    )

    if len(clusters_3) > 0:
        nx.draw_networkx_nodes(
            G,
            pos=pos,
            nodelist=list(clusters_3.index),
            node_color=color_sht,
            alpha=0.6,
            ax=ax,
            node_size=node_w.loc[clusters_3.index],
        )
        if node_labels == True:
            for node in G.nodes():
                if node in list(clusters_3.index):
                    labels[node] = node
                    labels_pos[node] = pos[node] + np.array([label_h, label_v])

    if node_labels == True:
        nx.draw_networkx_labels(
            G, pos=labels_pos, labels=labels, ax=ax, bbox=label_options, **font_options
        )

    ax.set_yticks([])
    ax.set_yticklabels([])
    for i in {"right", "left", "top", "bottom"}:
        side = ax.spines[i]
        side.set_visible(False)

    if title == "":
        if linkage == "DBHT":
            title = (
                "Planar Maximally Filtered Graph Asset Allocation("
                + codependence.capitalize()
                + " & "
                + kind
                + " layout)"
            )
        else:
            title = (
                "Minimun Spanning Tree Asset Allocation ("
                + codependence.capitalize()
                + " & "
                + kind
                + " layout)"
            )

    ax.set_title(title)

    try:
        fig.set_layout_engine(layout="constrained")
    except:
        pass

    return ax


def plot_clusters_network(
    returns,
    custom_cov=None,
    codependence="pearson",
    linkage="ward",
    opt_k_method="twodiff",
    k=None,
    max_k=10,
    bins_info="KN",
    alpha_tail=0.05,
    gs_threshold=0.5,
    leaf_order=True,
    seed=0,
    node_labels=True,
    node_size=2000,
    node_alpha=0.7,
    scale=10,
    subscale=5,
    font_size=10,
    title="",
    height=8,
    width=10,
    ax=None,
):
    r"""
    Create a network plot that show each cluster obtained from the dendrogram as an independent graph.

    Parameters
    ----------
    returns : DataFrame of shape (n_samples, n_assets)
        Assets returns DataFrame, where n_samples is the number of
        observations and n_assets is the number of assets.
    custom_cov : DataFrame or None, optional
        Custom covariance matrix, used when codependence parameter has value
        'custom_cov'. The default is None.
    codependence : str, can be {'pearson', 'spearman', 'abs_pearson', 'abs_spearman', 'distance', 'mutual_info', 'tail' or 'custom_cov'}
        The codependence or similarity matrix used to build the distance
        metric and clusters. The default is 'pearson'. Possible values are:

        - 'pearson': pearson correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{0.5(1-\rho^{pearson}_{i,j})}`.
        - 'spearman': spearman correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{0.5(1-\rho^{spearman}_{i,j})}`.
        - 'kendall': kendall correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{0.5(1-\rho^{kendall}_{i,j})}`.
        - 'gerber1': Gerber statistic 1 correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{0.5(1-\rho^{gerber1}_{i,j})}`.
        - 'gerber2': Gerber statistic 2 correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{0.5(1-\rho^{gerber2}_{i,j})}`.
        - 'abs_pearson': absolute value pearson correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{(1-|\rho_{i,j}|)}`.
        - 'abs_spearman': absolute value spearman correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{(1-|\rho_{i,j}|)}`.
        - 'abs_kendall': absolute value kendall correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{(1-|\rho^{kendall}_{i,j}|)}`.
        - 'distance': distance correlation matrix. Distance formula :math:`D_{i,j} = \sqrt{(1-\rho^{distance}_{i,j})}`.
        - 'mutual_info': mutual information matrix. Distance used is variation information matrix.
        - 'tail': lower tail dependence index matrix. Dissimilarity formula :math:`D_{i,j} = -\log{\lambda_{i,j}}`.
        - 'custom_cov': use custom correlation matrix based on the custom_cov parameter. Distance formula: :math:`D_{i,j} = \sqrt{0.5(1-\rho^{pearson}_{i,j})}`.

    linkage : string, optional
        Linkage method of hierarchical clustering, see `linkage <https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html?highlight=linkage#scipy.cluster.hierarchy.linkage>`_ for more details.
        The default is 'ward'. Possible values are:

        - 'single'.
        - 'complete'.
        - 'average'.
        - 'weighted'.
        - 'centroid'.
        - 'median'.
        - 'ward'.
        - 'DBHT': Direct Bubble Hierarchical Tree.

    opt_k_method : str
        Method used to calculate the optimum number of clusters.
        The default is 'twodiff'. Possible values are:

        - 'twodiff': two difference gap statistic.
        - 'stdsil': standarized silhouette score.

    k : int, optional
        Number of clusters. This value is took instead of the optimal number
        of clusters calculated with the two difference gap statistic.
        The default is None.
    max_k : int, optional
        Max number of clusters used by the two difference gap statistic
        to find the optimal number of clusters. The default is 10.
    bins_info: int or str
        Number of bins used to calculate variation of information. The default
        value is 'KN'. Possible values are:

        - 'KN': Knuth's choice method. See more in `knuth_bin_width <https://docs.astropy.org/en/stable/api/astropy.stats.knuth_bin_width.html>`_.
        - 'FD': Freedman–Diaconis' choice method. See more in `freedman_bin_width <https://docs.astropy.org/en/stable/api/astropy.stats.freedman_bin_width.html>`_.
        - 'SC': Scotts' choice method. See more in `scott_bin_width <https://docs.astropy.org/en/stable/api/astropy.stats.scott_bin_width.html>`_.
        - 'HGR': Hacine-Gharbi and Ravier' choice method.
        - int: integer value chosen by user.

    alpha_tail : float, optional
        Significance level for lower tail dependence index. The default is 0.05.
    gs_threshold : float, optional
        Gerber statistic threshold. The default is 0.5.
    leaf_order : bool, optional
        Indicates if the cluster are ordered so that the distance between
        successive leaves is minimal. The default is True.
    seed : int, optional
        Seed for networkx spring layout. The default value is 0.
    node_labels : bool, optional
        Specify if node lables are visible. The default value is True.
    node_size : float, optional
        Size of the node. The default value is 2000.
    node_alpha : float, optional
        Alpha parameter or transparency of nodes. The default value is 0.7.
    scale : float, optional
        Scale of whole graph. The default value is 10.
    subscale : float, optional
        Scale of clusters. The default value is 5.
    font_size : float, optional
        Font size of node labels. The default value is 12.
    title : str, optional
        Title of the chart. The default is "".
    height : float, optional
        Height of the image in inches. The default is 8.
    width : float, optional
        Width of the image in inches. The default is 10.
    ax : matplotlib axis, optional
        If provided, plot on this axis. The default is None.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    Returns
    -------
    ax : matplotlib axis
        Returns the Axes object with the plot for further tweaking.

    Example
    -------
    ::

        ax = rp.plot_clusters_network(returns=Y,
                                      codependence="pearson",
                                      linkage="ward",
                                      k=None,
                                      max_k=10,
                                      ax=None)

    .. image:: images/Assets_Clusters_Network.png


    """
    if not isinstance(returns, pd.DataFrame):
        raise ValueError("returns must be a DataFrame")

    if custom_cov is not None:
        if isinstance(custom_cov, pd.DataFrame):
            if custom_cov.index.tolist() != custom_cov.columns.tolist():
                raise ValueError(
                    "custom_cov must be a square DataFrame with samen labels in indexes and columns."
                )
            elif returns.columns.tolist() != custom_cov.index.tolist():
                raise ValueError("returns and custom_cov must have the same assets.")
        else:
            raise ValueError("custom_cov must be a square DataFrame.")

    if ax is None:
        fig = plt.gcf()
        ax = fig.gca()
        fig.set_figwidth(width)
        fig.set_figheight(height)
    else:
        fig = ax.get_figure()

    labels = np.array(returns.columns.tolist())

    # Calculating codependence matrix and distance metric
    codep, dist = af.codep_dist(
        returns=returns,
        custom_cov=custom_cov,
        codependence=codependence,
        bins_info=bins_info,
        alpha_tail=alpha_tail,
        gs_threshold=gs_threshold,
    )
    # Hierarchical clustering
    dist = dist.to_numpy()
    dist = pd.DataFrame(dist, columns=codep.columns, index=codep.index)
    if linkage == "DBHT":
        # different choices for D, S give different outputs!
        D = dist.to_numpy()  # dissimilarity matrix
        if codependence in {"pearson", "spearman", "custom_cov"}:
            S = (1 - dist**2).to_numpy()
        else:
            S = codep.copy().to_numpy()  # similarity matrix
        (_, Rpm, _, _, _, clustering) = db.DBHTs(
            D, S, leaf_order=leaf_order
        )  # DBHT clustering
    else:
        p_dist = squareform(dist, checks=False)
        clustering = hr.linkage(p_dist, method=linkage, optimal_ordering=leaf_order)

    # optimal number of clusters
    if k is None:
        if opt_k_method == "twodiff":
            k, clustering_inds = af.two_diff_gap_stat(dist, clustering, max_k)
        elif opt_k_method == "stdsil":
            k, clustering_inds = af.std_silhouette_score(dist, clustering, max_k)
        else:
            raise ValueError("The only opt_k_method available are twodiff and stdsil")
    else:
        clustering_inds = hr.fcluster(clustering, k, criterion="maxclust")

    clusters = {i: [] for i in np.unique(clustering_inds)}
    for i, v in enumerate(clustering_inds):
        clusters[v].append(labels[i])

    # Calculating adjacency matrix based on hierarchical clustering
    D = ct.clusters_matrix(
        returns,
        codependence=codependence,
        linkage=linkage,
        k=k,
        max_k=max_k,
        bins_info=bins_info,
        alpha_tail=alpha_tail,
        gs_threshold=gs_threshold,
        leaf_order=leaf_order,
    )

    # Build Adjacency matrix
    MAdj = pd.DataFrame(D, index=labels, columns=labels)
    G = nx.from_pandas_adjacency(MAdj)

    # Build clusters super nodes positions
    supergraph = nx.cycle_graph(len(clusters))
    superpos = nx.spring_layout(supergraph, scale=scale, seed=int(seed))

    # Use the "supernode" positions as the center of each node cluster
    centers = list(superpos.values())
    pos = {}
    for center, key in zip(centers, clusters.keys()):
        pos.update(
            nx.spring_layout(
                nx.subgraph(G, clusters[key]),
                center=center,
                scale=subscale,
                seed=int(seed),
            )
        )

    colors = af.color_list(k)

    label_options = {"ec": "k", "fc": "white", "alpha": 0.7}

    node_options = {
        "node_size": node_size,
        "alpha": node_alpha,
    }
    font_options = {
        "font_size": font_size,
        "font_color": "k",
    }
    # Nodes colored by cluster
    for key, clr in zip(clusters.keys(), colors):
        nx.draw_networkx_nodes(
            G, pos=pos, nodelist=clusters[key], node_color=clr, **node_options
        )
    nx.draw_networkx_edges(G, pos=pos, ax=ax, edge_color="grey")

    if node_labels == True:
        nx.draw_networkx_labels(G, pos=pos, ax=ax, bbox=label_options, **font_options)

    ax.set_yticks([])
    ax.set_yticklabels([])
    for i in {"right", "left", "top", "bottom"}:
        side = ax.spines[i]
        side.set_visible(False)

    if title == "":
        title = (
            "Cluster Network ("
            + codependence.capitalize()
            + " & "
            + linkage
            + " linkage)"
        )

    ax.set_title(title)

    try:
        fig.set_layout_engine(layout="constrained")
    except:
        pass

    return ax


def plot_clusters_network_allocation(
    returns,
    w,
    custom_cov=None,
    codependence="pearson",
    linkage="ward",
    opt_k_method="twodiff",
    k=None,
    max_k=10,
    bins_info="KN",
    alpha_tail=0.05,
    gs_threshold=0.5,
    leaf_order=True,
    seed=0,
    node_labels=True,
    max_node_size=2000,
    color_lng="tab:blue",
    color_sht="tab:red",
    scale=10,
    subscale=5,
    label_v=1.5,
    label_h=0,
    font_size=10,
    title="",
    height=8,
    width=10,
    ax=None,
):
    r"""
    Create a network plot that show each cluster obtained from the dendrogram as
    an independent graph. The size of the nodes and color represents the amount invested and direction
    (long-short) respectively.

    Parameters
    ----------
    returns : DataFrame of shape (n_samples, n_assets)
        Assets returns DataFrame, where n_samples is the number of
        observations and n_assets is the number of assets.
    w : DataFrame or Series of shape (n_assets, 1)
        Portfolio weights, where n_assets is the number of assets.
    custom_cov : DataFrame or None, optional
        Custom covariance matrix, used when codependence parameter has value
        'custom_cov'. The default is None.
    codependence : str, can be {'pearson', 'spearman', 'abs_pearson', 'abs_spearman', 'distance', 'mutual_info', 'tail' or 'custom_cov'}
        The codependence or similarity matrix used to build the distance
        metric and clusters. The default is 'pearson'. Possible values are:

        - 'pearson': pearson correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{0.5(1-\rho^{pearson}_{i,j})}`.
        - 'spearman': spearman correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{0.5(1-\rho^{spearman}_{i,j})}`.
        - 'kendall': kendall correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{0.5(1-\rho^{kendall}_{i,j})}`.
        - 'gerber1': Gerber statistic 1 correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{0.5(1-\rho^{gerber1}_{i,j})}`.
        - 'gerber2': Gerber statistic 2 correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{0.5(1-\rho^{gerber2}_{i,j})}`.
        - 'abs_pearson': absolute value pearson correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{(1-|\rho_{i,j}|)}`.
        - 'abs_spearman': absolute value spearman correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{(1-|\rho_{i,j}|)}`.
        - 'abs_kendall': absolute value kendall correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{(1-|\rho^{kendall}_{i,j}|)}`.
        - 'distance': distance correlation matrix. Distance formula :math:`D_{i,j} = \sqrt{(1-\rho^{distance}_{i,j})}`.
        - 'mutual_info': mutual information matrix. Distance used is variation information matrix.
        - 'tail': lower tail dependence index matrix. Dissimilarity formula :math:`D_{i,j} = -\log{\lambda_{i,j}}`.
        - 'custom_cov': use custom correlation matrix based on the custom_cov parameter. Distance formula: :math:`D_{i,j} = \sqrt{0.5(1-\rho^{pearson}_{i,j})}`.

    linkage : string, optional
        Linkage method of hierarchical clustering, see `linkage <https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html?highlight=linkage#scipy.cluster.hierarchy.linkage>`_ for more details.
        The default is 'ward'. Possible values are:

        - 'single'.
        - 'complete'.
        - 'average'.
        - 'weighted'.
        - 'centroid'.
        - 'median'.
        - 'ward'.
        - 'DBHT': Direct Bubble Hierarchical Tree.

    bins_info: int or str
        Number of bins used to calculate variation of information. The default
        value is 'KN'. Possible values are:

        - 'KN': Knuth's choice method. See more in `knuth_bin_width <https://docs.astropy.org/en/stable/api/astropy.stats.knuth_bin_width.html>`_.
        - 'FD': Freedman–Diaconis' choice method. See more in `freedman_bin_width <https://docs.astropy.org/en/stable/api/astropy.stats.freedman_bin_width.html>`_.
        - 'SC': Scotts' choice method. See more in `scott_bin_width <https://docs.astropy.org/en/stable/api/astropy.stats.scott_bin_width.html>`_.
        - 'HGR': Hacine-Gharbi and Ravier' choice method.
        - int: integer value chosen by user.

    alpha_tail : float, optional
        Significance level for lower tail dependence index. The default is 0.05.
    gs_threshold : float, optional
        Gerber statistic threshold. The default is 0.5.
    leaf_order : bool, optional
        Indicates if the cluster are ordered so that the distance between
        successive leaves is minimal. The default is True.
    seed : int, optional
        Seed for networkx spring layout. The default value is 0.
    node_labels : bool, optional
        Specify if node lables are visible. The default value is True.
    max_node_size : float, optional
        Size of the node with maximum weight in absolute value. The default
        value is 2000.
    color_lng : str, optional
        Color of assets with long positions. The default value is 'tab:blue'.
    color_sht : str, optional
        Color of assets with short positions. The default value is 'tab:red'.
    scale : float, optional
        Scale of whole graph. The default value is 10.
    subscale : float, optional
        Scale of clusters. The default value is 5.
    label_v : float, optional
        Vertical distance the label is offset from the center. The default value is 1.5.
    label_h : float, optional
        Horizontal distance the label is offset from the center. The default value is 0.
    font_size : float, optional
        Font size of node labels. The default value is 10.
    title : str, optional
        Title of the chart. The default is "".
    height : float, optional
        Height of the image in inches. The default is 8.
    width : float, optional
        Width of the image in inches. The default is 10.
    ax : matplotlib axis, optional
        If provided, plot on this axis. The default is None.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    Returns
    -------
    ax : matplotlib axis
        Returns the Axes object with the plot for further tweaking.

    Example
    -------
    ::

        ax = rp.plot_clusters_network_allocation(returns=Y,
                                                 w=w1,
                                                 codependence="pearson",
                                                 linkage="ward",
                                                 k=None,
                                                 max_k=10,
                                                 ax=None)

    .. image:: images/Assets_Clusters_Network_Allocation.png


    """
    if not isinstance(returns, pd.DataFrame):
        raise ValueError("returns must be a DataFrame")

    if not isinstance(w, pd.DataFrame):
        if isinstance(w, pd.Series):
            w_ = w.to_frame()
        else:
            raise ValueError("w must be a one column DataFrame or Series")
    else:
        if w.shape[0] == 1:
            w_ = w.T.copy()
        elif w.shape[1] == 1:
            w_ = w.copy()
        else:
            raise ValueError("w must be a one column DataFrame or Series")

    if returns.columns.tolist() != w_.index.tolist():
        raise ValueError("returns and w must have same assets.")

    if custom_cov is not None:
        if isinstance(custom_cov, pd.DataFrame):
            if custom_cov.index.tolist() != custom_cov.columns.tolist():
                raise ValueError(
                    "custom_cov must be a square DataFrame with samen labels in indexes and columns."
                )
            elif returns.columns.tolist() != custom_cov.index.tolist():
                raise ValueError("returns and custom_cov must have the same assets.")
        else:
            raise ValueError("custom_cov must be a square DataFrame.")

    if ax is None:
        fig = plt.gcf()
        ax = fig.gca()
        fig.set_figwidth(width)
        fig.set_figheight(height)
    else:
        fig = ax.get_figure()

    labels = np.array(returns.columns.tolist())

    # Calculating codependence matrix and distance metric
    codep, dist = af.codep_dist(
        returns=returns,
        custom_cov=custom_cov,
        codependence=codependence,
        bins_info=bins_info,
        alpha_tail=alpha_tail,
        gs_threshold=gs_threshold,
    )

    # Hierarchical clustering
    dist = dist.to_numpy()
    dist = pd.DataFrame(dist, columns=codep.columns, index=codep.index)
    if linkage == "DBHT":
        # different choices for D, S give different outputs!
        D = dist.to_numpy()  # dissimilarity matrix
        if codependence in {"pearson", "spearman", "custom_cov"}:
            S = (1 - dist**2).to_numpy()
        else:
            S = codep.copy().to_numpy()  # similarity matrix
        (_, Rpm, _, _, _, clustering) = db.DBHTs(
            D, S, leaf_order=leaf_order
        )  # DBHT clustering
    else:
        p_dist = squareform(dist, checks=False)
        clustering = hr.linkage(p_dist, method=linkage, optimal_ordering=leaf_order)

    # optimal number of clusters
    if k is None:
        if opt_k_method == "twodiff":
            k, clustering_inds = af.two_diff_gap_stat(dist, clustering, max_k)
        elif opt_k_method == "stdsil":
            k, clustering_inds = af.std_silhouette_score(dist, clustering, max_k)
        else:
            raise ValueError("The only opt_k_method available are twodiff and stdsil")
    else:
        clustering_inds = hr.fcluster(clustering, k, criterion="maxclust")

    clusters = {i: [] for i in np.unique(clustering_inds)}
    for i, v in enumerate(clustering_inds):
        clusters[v].append(labels[i])

    # Calculating adjacency matrix based on hierarchical clustering
    D = ct.clusters_matrix(
        returns,
        codependence=codependence,
        linkage=linkage,
        k=k,
        max_k=max_k,
        bins_info=bins_info,
        alpha_tail=alpha_tail,
        gs_threshold=gs_threshold,
        leaf_order=leaf_order,
    )

    # Build Adjacency matrix
    MAdj = pd.DataFrame(D, index=labels, columns=labels)
    G = nx.from_pandas_adjacency(MAdj)

    # Build clusters super nodes positions
    supergraph = nx.cycle_graph(len(clusters))
    superpos = nx.spring_layout(supergraph, scale=scale, seed=int(seed))

    # Use the "supernode" positions as the center of each node cluster
    centers = list(superpos.values())
    pos = {}
    for center, key in zip(centers, clusters.keys()):
        pos.update(
            nx.spring_layout(
                nx.subgraph(G, clusters[key]),
                center=center,
                scale=subscale,
                seed=int(seed),
            )
        )

    label_options = {"ec": "k", "fc": "white", "alpha": 0.7}

    font_options = {
        "font_size": font_size,
        "font_color": "k",
    }

    clusters_1 = set(w_.loc[w_.iloc[:, 0] > 1e-6].index)
    clusters_2 = set(w_.loc[(w_.iloc[:, 0] <= 1e-6) & (w_.iloc[:, 0] >= -1e-6)].index)
    clusters_3 = set(w_.loc[w_.iloc[:, 0] < -1e-6].index)
    node_w = np.abs(w_) * max_node_size / np.abs(w_).sum().item()

    # Nodes colored by cluster
    for key in clusters.keys():
        nodes_pos = list(clusters_1.intersection(clusters[key]))
        nodes_0 = list(clusters_2.intersection(clusters[key]))
        nx.draw_networkx_nodes(
            G,
            pos=pos,
            nodelist=nodes_pos,
            node_color=color_lng,
            alpha=0.8,
            ax=ax,
            node_size=node_w.loc[nodes_pos],
        )
        nx.draw_networkx_nodes(
            G,
            pos=pos,
            nodelist=nodes_0,
            node_color="lightgrey",
            alpha=0.5,
            ax=ax,
            node_size=100,
        )
        if len(clusters_3) > 0:
            nodes_neg = list(clusters_3.intersection(clusters[key]))
            nx.draw_networkx_nodes(
                G,
                pos=pos,
                nodelist=nodes_neg,
                node_color=color_sht,
                alpha=0.6,
                ax=ax,
                node_size=node_w.loc[nodes_neg],
            )

    nx.draw_networkx_edges(G, pos=pos, ax=ax, edge_color="grey")

    if node_labels == True:
        label_w = list(w_.loc[np.abs(w_.iloc[:, 0]) >= 1e-6].index)
        labels = {}
        labels_pos = {}
        for node in G.nodes():
            if node in label_w:
                labels[node] = node
                labels_pos[node] = pos[node] + np.array([label_h, label_v])

    if node_labels == True:
        nx.draw_networkx_labels(
            G, pos=labels_pos, labels=labels, ax=ax, bbox=label_options, **font_options
        )

    ax.set_yticks([])
    ax.set_yticklabels([])
    for i in {"right", "left", "top", "bottom"}:
        side = ax.spines[i]
        side.set_visible(False)

    if title == "":
        title = (
            "Clusters Network Asset Allocation ("
            + codependence.capitalize()
            + " & "
            + linkage
            + " linkage)"
        )

    ax.set_title(title)

    try:
        fig.set_layout_engine(layout="constrained")
    except:
        pass

    return ax


def plot_BrinsonAttribution(
    prices,
    w,
    wb,
    start,
    end,
    asset_classes,
    classes_col,
    method="nearest",
    sector="Total",
    height=6,
    width=10,
    ax=None,
):
    r"""
    Creates a plot with the Brinson Performance Attribution specified by the
    sector parameter.

    Parameters
    ----------
    prices : DataFrame of shape (n_samples, n_assets)
        Assets prices DataFrame, where n_samples is the number of
        observations and n_assets is the number of assets.
    w : DataFrame or Series of shape (n_assets, 1)
        A portfolio specified by the user.
    wb : DataFrame or Series of shape (n_assets, 1)
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
        See `get_indexer <https://pandas.pydata.org/docs/reference/api/pandas.Index.get_indexer.html#pandas.Index.get_indexer>`_ for more details.
    sector : str
        Is the sector or class for which the function will plot the Brinson
        performance attribution. Possible values are all classes in the set
        of classes specified by classes_col parameter and 'Total' for aggregate
        performance attribution. Default value is 'Total'.
    height : float, optional
        Height of the image in inches. The default is 8.
    width : float, optional
        Width of the image in inches. The default is 10.
    ax : matplotlib axis, optional
        If provided, plot on this axis. The default is None.

    Raises
    ------
    ValueError
        When the value cannot be calculated.

    Returns
    -------
    ax :  matplotlib axis.
        Returns the Axes object with the plot for further tweaking.

    Example
    -------
    ::

        ax = plot_BrinsonAttribution(
            prices=data,
            w=w,
            wb=wb,
            start='2019-01-07',
            end='2019-12-06',
            asset_classes=asset_classes,
            classes_col='Industry',
            method='nearest',
            sector='Total',
            height=6,
            width=10,
            ax=None
            )

    .. image:: images/BrinAttr_plot.png


    """

    if ax is None:
        fig = plt.gcf()
        ax = fig.gca()
        fig.set_figwidth(width)
        fig.set_figheight(height)
    else:
        fig = ax.get_figure()

    BrinAttr, (start_, end_) = rk.BrinsonAttribution(
        prices=prices,
        w=w,
        wb=wb,
        start=start,
        end=end,
        asset_classes=asset_classes,
        classes_col=classes_col,
        method=method,
    )

    if sector not in BrinAttr.columns.tolist():
        raise ValueError("Sector is not in asset_classes or it is not the Total")
    else:
        labels = [
            "Asset Allocation",
            "Security Selection",
            "Interaction",
            "Total Excess Return",
        ]

        ax.barh(
            3,
            BrinAttr.loc[labels[0], sector],
            align="center",
            color="fuchsia",
            edgecolor="black",
            alpha=0.3,
        )
        ax.barh(
            2,
            BrinAttr.loc[labels[1], sector],
            align="center",
            color="orange",
            edgecolor="black",
            alpha=0.3,
        )
        ax.barh(
            1,
            BrinAttr.loc[labels[2], sector],
            align="center",
            color="lime",
            edgecolor="black",
            alpha=0.3,
        )
        if BrinAttr.loc[labels[3], sector] < 0:
            ax.barh(
                0,
                BrinAttr.loc[labels[3], sector],
                align="center",
                color="r",
                edgecolor="black",
                alpha=0.3,
            )
        else:
            ax.barh(
                0,
                BrinAttr.loc[labels[3], sector],
                align="center",
                color="b",
                edgecolor="black",
                alpha=0.3,
            )
        ax.set_yticks([3, 2, 1, 0])
        ax.set_yticklabels(labels)
        ax.set_xlabel("Excess Return")
        ax.set_title(
            sector + " Performance Attribution Chart (" + start_ + " to " + end_ + ")"
        )

        ax.set_xticklabels(["{:3.2%}".format(x) for x in ax.get_xticks()])

        r = plt.gcf().canvas.get_renderer()
        transf = ax.transData.inverted()
        sizes = [
            BrinAttr.loc[labels[3], sector],
            BrinAttr.loc[labels[2], sector],
            BrinAttr.loc[labels[1], sector],
            BrinAttr.loc[labels[0], sector],
        ]
        sizes2 = ["{:.2%}".format(x) for x in sizes]

        for i, v in enumerate(sizes):
            t = ax.text(v, i, sizes2[i], color="black")
            bb = t.get_window_extent(renderer=r)
            bb = bb.transformed(transf)
            w_text = bb.width
            x_text = bb.x0
            y_text = bb.y0
            if v >= 0:
                t.set_position((x_text + w_text * 0.2, y_text))
            else:
                t.set_position((x_text - w_text * 1.2, y_text))

        ax.set_xlim(
            min(sizes) - (max(sizes) - min(sizes)) / 7,
            max(sizes) + (max(sizes) - min(sizes)) / 7,
        )

        ax.grid(linestyle=":")

        try:
            fig.set_layout_engine(layout="constrained")
        except:
            pass

    return ax
