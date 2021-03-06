import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib import cm
import scipy.stats as st
import riskfolio.RiskFunctions as rk

__all__ = [
    "plot_series",
    "plot_frontier",
    "plot_pie",
    "plot_frontier_area",
    "plot_risk_con",
    "plot_hist",
    "plot_drawdown",
    "plot_table",
]

rm_names = [
    "Standard Deviation",
    "Mean Absolute Deviation",
    "Semi Standard Deviation",
    "Value at Risk",
    "Conditional Value at Risk",
    "Entropic Value at Risk",
    "Worst Realization",
    "First Lower Partial Moment",
    "Second Lower Partial Moment",
    "Max Drawdown",
    "Average Drawdown",
    "Drawdown at Risk",
    "Conditional Drawdown at Risk",
    "Entropic Drawdown at Risk",
    "Ulcer Index",
]

rmeasures = [
    "MV",
    "MAD",
    "MSV",
    "VaR",
    "CVaR",
    "EVaR",
    "WR",
    "FLPM",
    "SLPM",
    "MDD",
    "ADD",
    "DaR",
    "CDaR",
    "EDaR",
    "UCI",
]


def plot_series(returns, w, cmap="tab20", height=6, width=10, ax=None):
    r"""
    Create a chart with the compound cumulated of the portfolios.

    Parameters
    ----------
    returns : DataFrame
        Assets returns.
    w : DataFrame of shape (n_assets, n_portfolios)
        Portfolio weights.
    cmap : cmap, optional
        Colorscale, represente the risk adjusted return ratio.
        The default is 'tab20'.
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

        ax = plf.plot_series(returns=Y, w=ws, cmap='tab20', height=6, width=10, ax=None)

    .. image:: images/Port_Series.png


    """
    if not isinstance(returns, pd.DataFrame):
        raise ValueError("returns must be a DataFrame")

    if not isinstance(w, pd.DataFrame):
        raise ValueError("w must be a DataFrame")

    if returns.shape[1] != w.shape[0]:
        a1 = str(returns.shape)
        a2 = str(w.shape)
        raise ValueError("shapes " + a1 + " and " + a2 + " not aligned")

    if ax is None:
        ax = plt.gca()
        fig = plt.gcf()
        fig.set_figwidth(width)
        fig.set_figheight(height)

    ax.grid(linestyle=":")
    title = "Historical Compounded Cumulative Returns"
    ax.set_title(title)

    labels = w.columns.tolist()

    colormap = cm.get_cmap(cmap)
    colormap = colormap(np.linspace(0, 1, 20))

    if cmap == "gist_rainbow":
        colormap = colormap[::-1]

    cycle = plt.cycler("color", colormap)
    ax.set_prop_cycle(cycle)

    X = w.columns.tolist()
    index = returns.index.tolist()

    for i in range(len(X)):
        a = np.array(returns, ndmin=2) @ np.array(w[X[i]], ndmin=2).T
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

    fig = plt.gcf()
    fig.tight_layout()

    return ax


def plot_frontier(
    w_frontier,
    mu,
    cov=None,
    returns=None,
    rm="MV",
    rf=0,
    alpha=0.05,
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
    cov : DataFrame of shape (n_features, n_features)
        Covariance matrix, where n_features is the number of features.
    returns : DataFrame of shape (n_samples, n_features)
        Features matrix, where n_samples is the number of samples and
        n_features is the number of features.
    rm : str, optional
        The risk measure used to estimate the frontier.
        The default is 'MV'. Posible values are:

        - 'MV': Standard Deviation.
        - 'MAD': Mean Absolute Deviation.
        - 'MSV': Semi Standard Deviation.
        - 'FLPM': First Lower Partial Moment (Omega Ratio).
        - 'SLPM': Second Lower Partial Moment (Sortino Ratio).
        - 'CVaR': Conditional Value at Risk.
        - 'EVaR': Conditional Value at Risk.
        - 'WR': Worst Realization (Minimax)
        - 'MDD': Maximum Drawdown of uncompounded returns (Calmar Ratio).
        - 'ADD': Average Drawdown of uncompounded cumulative returns.
        - 'DaR': Drawdown at Risk of uncompounded cumulative returns.
        - 'CDaR': Conditional Drawdown at Risk of uncompounded cumulative returns.
        - 'EDaR': Entropic Drawdown at Risk of uncompounded cumulative returns.
        - 'UCI': Ulcer Index of uncompounded cumulative returns.

    rf : float, optional
        Risk free rate or minimum aceptable return. The default is 0.
    alpha : float, optional
        Significante level of VaR, CVaR, EVaR, DaR and CDaR.
        The default is 0.05.
    cmap : cmap, optional
        Colorscale, represente the risk adjusted return ratio.
        The default is 'viridis'.
    w : DataFrame of shape (n_assets, 1), optional
        A portfolio specified by the user. The default is None.
    label : str, optional
        Name of portfolio that appear on plot legend.
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

        ax = plf.plot_frontier(w_frontier=ws, mu=mu, cov=cov, returns=returns,
                               rm=rm, rf=0, alpha=0.05, cmap='viridis', w=w1,
                               label=label, marker='*', s=16, c='r',
                               height=6, width=10, t_factor=252, ax=None)

    .. image:: images/MSV_Frontier.png


    """

    if not isinstance(w_frontier, pd.DataFrame):
        raise ValueError("w_frontier must be a DataFrame")

    if not isinstance(mu, pd.DataFrame):
        raise ValueError("mu must be a DataFrame")

    if not isinstance(cov, pd.DataFrame):
        raise ValueError("cov must be a DataFrame")

    if not isinstance(returns, pd.DataFrame):
        raise ValueError("returns must be a DataFrame")

    if returns.shape[1] != w_frontier.shape[0]:
        a1 = str(returns.shape)
        a2 = str(w_frontier.shape)
        raise ValueError("shapes " + a1 + " and " + a2 + " not aligned")

    if w is not None:
        if not isinstance(w, pd.DataFrame):
            raise ValueError("w must be a DataFrame")

        if w.shape[1] > 1 and w.shape[0] == 0:
            w = w.T
        elif w.shape[1] > 1 and w.shape[0] > 0:
            raise ValueError("w must be a column DataFrame")

        if returns.shape[1] != w.shape[0]:
            a1 = str(returns.shape)
            a2 = str(w.shape)
            raise ValueError("shapes " + a1 + " and " + a2 + " not aligned")

    if ax is None:
        ax = plt.gca()
        fig = plt.gcf()
        fig.set_figwidth(width)
        fig.set_figheight(height)

    mu_ = np.array(mu, ndmin=2)

    ax.set_ylabel("Expected Return")
    item = rmeasures.index(rm)
    x_label = rm_names[item] + " (" + rm + ")"
    ax.set_xlabel("Expected Risk - " + x_label)

    title = "Efficient Frontier Mean - " + x_label
    ax.set_title(title)

    X1 = []
    Y1 = []
    Z1 = []

    for i in range(w_frontier.shape[1]):
        weights = np.array(w_frontier.iloc[:, i], ndmin=2).T
        risk = rk.Sharpe_Risk(
            weights, cov=cov, returns=returns, rm=rm, rf=rf, alpha=alpha
        )

        ret = mu_ @ weights
        ret = ret.item() * t_factor

        if rm not in ["MDD", "ADD", "CDaR", "EDaR", "UCI"]:
            risk = risk * t_factor ** 0.5

        ratio = (ret - rf) / risk

        X1.append(risk)
        Y1.append(ret)
        Z1.append(ratio)

    ax1 = ax.scatter(X1, Y1, c=Z1, cmap=cmap)

    if w is not None:
        X2 = []
        Y2 = []
        for i in range(w.shape[1]):
            weights = np.array(w.iloc[:, i], ndmin=2).T
            risk = rk.Sharpe_Risk(
                weights, cov=cov, returns=returns, rm=rm, rf=rf, alpha=alpha
            )
            ret = mu_ @ weights
            ret = ret.item() * t_factor

            if rm not in ["MDD", "ADD", "CDaR", "EDaR", "UCI"]:
                risk = risk * t_factor ** 0.5

            ratio = (ret - rf) / risk

            X2.append(risk)
            Y2.append(ret)

        ax.scatter(X2, Y2, marker=marker, s=s ** 2, c=c, label=label)
        ax.legend(loc="upper left")

    xmin = np.min(X1) - np.abs(np.max(X1) - np.min(X1)) * 0.1
    xmax = np.max(X1) + np.abs(np.max(X1) - np.min(X1)) * 0.1
    ymin = np.min(Y1) - np.abs(np.max(Y1) - np.min(Y1)) * 0.1
    ymax = np.max(Y1) + np.abs(np.max(Y1) - np.min(Y1)) * 0.1

    ax.set_ylim(ymin, ymax)
    ax.set_xlim(xmin, xmax)

    ax.xaxis.set_major_locator(plt.AutoLocator())

    ticks_loc = ax.get_yticks().tolist()
    ax.set_yticks(ax.get_yticks().tolist())
    ax.set_yticklabels(["{:.2%}".format(x) for x in ticks_loc])
    ticks_loc = ax.get_xticks().tolist()
    ax.set_xticks(ax.get_xticks().tolist())
    ax.set_xticklabels(["{:.2%}".format(x) for x in ticks_loc])

    ax.tick_params(axis="y", direction="in")
    ax.tick_params(axis="x", direction="in")

    ax.grid(linestyle=":")

    colorbar = ax.figure.colorbar(ax1)
    colorbar.set_label("Risk Adjusted Return Ratio")

    fig = plt.gcf()
    fig.tight_layout()

    return ax


def plot_pie(
    w, title="", others=0.05, nrow=25, cmap="tab20", height=6, width=8, ax=None
):
    r"""
    Create a pie chart with portfolio weights.

    Parameters
    ----------
    w : DataFrame of shape (n_assets, 1)
        Portfolio weights.
    title : str, optional
        Title of the chart. The default is ''.
    others : float, optional
        Percentage of others section. The default is 0.05.
    nrow : int, optional
        Number of rows of the legend. The default is 25.
    cmap : cmap, optional
        Color scale used to plot each asset weight.
        The default is 'tab20'.
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

        ax = plf.plot_pie(w=w1, title='Portafolio', height=6, width=10,
                          cmap="tab20", ax=None)

    .. image:: images/Pie_Chart.png


    """

    if not isinstance(w, pd.DataFrame):
        raise ValueError("w must be a DataFrame")

    if w.shape[1] > 1 and w.shape[0] == 0:
        w = w.T
    elif w.shape[1] > 1 and w.shape[0] > 0:
        raise ValueError("w must be a column DataFrame")

    if ax is None:
        ax = plt.gca()
        fig = plt.gcf()
        fig.set_figwidth(width)
        fig.set_figheight(height)

    labels = w.index.tolist()
    sizes = w.iloc[:, 0].tolist()
    abs_sizes = [np.absolute(s) for s in sizes]
    sizes2 = pd.DataFrame([labels, abs_sizes, sizes]).T
    sizes2.columns = ["labels", "abs_values", "values"]
    sizes2 = sizes2.sort_values(by=["abs_values"], ascending=False)
    sizes2.index = [i for i in range(0, len(labels))]
    sizes3 = sizes2.cumsum()
    sizes3["abs_values"] = sizes3["abs_values"] / sizes3["abs_values"].max()
    l = sizes3[sizes3["abs_values"] >= 1 - others].index.tolist()[0]

    a1 = sizes2["abs_values"].sum() - sizes2[sizes2.index <= l]["abs_values"].sum()
    a2 = sizes2["values"].sum() - sizes2[sizes2.index <= l]["values"].sum()
    item = pd.DataFrame(["Others", a1, a2]).T
    item.columns = ["labels", "abs_values", "values"]
    sizes2 = sizes2[sizes2.index <= l]
    sizes2 = sizes2.append(item)

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
    colormap = colormap(np.linspace(0, 1, 20))

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

    fig = plt.gcf()
    fig.tight_layout()

    return ax


def plot_frontier_area(w_frontier, nrow=25, cmap="tab20", height=6, width=10, ax=None):
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

        ax = plf.plot_frontier_area(w_frontier=ws, cmap="tab20", height=6,
                                    width=10, ax=None)

    .. image:: images/Area_Frontier.png


    """

    if not isinstance(w_frontier, pd.DataFrame):
        raise ValueError("w must be a DataFrame")

    if ax is None:
        ax = plt.gca()
        fig = plt.gcf()
        fig.set_figwidth(width)
        fig.set_figheight(height)

    ax.set_title("Efficient Frontier's Assets Structure")
    labels = w_frontier.index.tolist()

    colormap = cm.get_cmap(cmap)
    colormap = colormap(np.linspace(0, 1, 20))

    if cmap == "gist_rainbow":
        colormap = colormap[::-1]

    cycle = plt.cycler("color", colormap)
    ax.set_prop_cycle(cycle)

    X = w_frontier.columns.tolist()

    ax.stackplot(X, w_frontier, labels=labels, alpha=0.7, edgecolor="black")

    ax.set_ylim(0, 1)
    ax.set_xlim(0, len(X) - 1)

    ticks_loc = ax.get_yticks().tolist()
    ax.set_yticks(ax.get_yticks().tolist())
    ax.set_yticklabels(["{:3.2%}".format(x) for x in ticks_loc])
    ax.grid(linestyle=":")

    n = int(np.ceil(len(labels) / nrow))

    ax.legend(labels, loc="center left", bbox_to_anchor=(1, 0.5), ncol=n)

    fig = plt.gcf()
    fig.tight_layout()

    return ax


def plot_risk_con(
    w,
    cov=None,
    returns=None,
    rm="MV",
    rf=0,
    alpha=0.05,
    color="tab:blue",
    height=6,
    width=10,
    t_factor=252,
    ax=None,
):
    r"""
    Create a chart with the risk contribution per asset of the portfolio.

    Parameters
    ----------
    w : DataFrame of shape (n_assets, 1)
        Portfolio weights.
    cov : DataFrame of shape (n_features, n_features)
        Covariance matrix, where n_features is the number of features.
    returns : DataFrame of shape (n_samples, n_features)
        Features matrix, where n_samples is the number of samples and
        n_features is the number of features.
    rm : str, optional
        Risk measure used to estimate risk contribution.
        The default is 'MV'. Posible values are:

        - 'MV': Standard Deviation.
        - 'MAD': Mean Absolute Deviation.
        - 'MSV': Semi Standard Deviation.
        - 'FLPM': First Lower Partial Moment (Omega Ratio).
        - 'SLPM': Second Lower Partial Moment (Sortino Ratio).
        - 'CVaR': Conditional Value at Risk.
        - 'EVaR': Conditional Value at Risk.
        - 'WR': Worst Realization (Minimax)
        - 'MDD': Maximum Drawdown of uncompounded cumulative returns (Calmar Ratio).
        - 'ADD': Average Drawdown of uncompounded cumulative returns.
        - 'DaR': Drawdown at Risk of uncompounded cumulative returns.
        - 'CDaR': Conditional Drawdown at Risk of uncompounded cumulative returns.
        - 'UCI': Ulcer Index of uncompounded cumulative returns.

    rf : float, optional
        Risk free rate or minimum aceptable return. The default is 0.
    alpha : float, optional
        Significante level of VaR, CVaR and CDaR. The default is 0.05.
    color : str, optional
        Color used to plot each asset risk contribution.
        The default is 'tab:blue'.
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

        ax = plf.plot_risk_con(w=w2, cov=cov, returns=returns, rm=rm,
                               rf=0, alpha=0.05, color="tab:blue", height=6,
                               width=10, t_factor=252, ax=None)

    .. image:: images/Risk_Con.png


    """

    if not isinstance(w, pd.DataFrame):
        raise ValueError("w must be a DataFrame")

    if ax is None:
        ax = plt.gca()
        fig = plt.gcf()
        fig.set_figwidth(width)
        fig.set_figheight(height)

    item = rmeasures.index(rm)
    title = "Risk (" + rm_names[item] + ") Contribution per Asset"
    ax.set_title(title)

    X = w.index.tolist()

    RC = rk.Risk_Contribution(w, cov=cov, returns=returns, rm=rm, rf=rf, alpha=alpha)

    if rm not in ["MDD", "ADD", "CDaR", "EDaR", "UCI"]:
        RC = RC * t_factor ** 0.5

    ax.bar(X, RC, alpha=0.7, color=color, edgecolor="black")

    ax.set_xlim(-0.5, len(X) - 0.5)

    ticks_loc = ax.get_yticks().tolist()
    ax.set_yticks(ax.get_yticks())
    ax.set_yticklabels(["{:3.4%}".format(x) for x in ticks_loc])
    ax.grid(linestyle=":")

    fig = plt.gcf()
    fig.tight_layout()

    return ax


def plot_hist(returns, w, alpha=0.05, bins=50, height=6, width=10, ax=None):
    r"""
    Create a histogram of portfolio returns with the risk measures.

    Parameters
    ----------
    returns : DataFrame
        Assets returns.
    w : DataFrame of shape (n_assets, 1)
        Portfolio weights.
    alpha : float, optional
        Significante level of VaR, CVaR and EVaR. The default is 0.05.
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

        ax = plf.plot_hist(returns=Y, w=w1, alpha=0.05, bins=50, height=6,
                           width=10, ax=None)

    .. image:: images/Histogram.png


    """

    if not isinstance(returns, pd.DataFrame):
        raise ValueError("returns must be a DataFrame")

    if not isinstance(w, pd.DataFrame):
        raise ValueError("w must be a DataFrame")

    if w.shape[1] > 1 and w.shape[0] == 0:
        w = w.T
    elif w.shape[1] > 1 and w.shape[0] > 0:
        raise ValueError("w must be a  DataFrame")

    if returns.shape[1] != w.shape[0]:
        a1 = str(returns.shape)
        a2 = str(w.shape)
        raise ValueError("shapes " + a1 + " and " + a2 + " not aligned")

    if ax is None:
        ax = plt.gca()
        fig = plt.gcf()
        fig.set_figwidth(width)
        fig.set_figheight(height)

    a = np.array(returns, ndmin=2) @ np.array(w, ndmin=2)
    ax.set_title("Portfolio Returns Histogram")
    n, bins1, patches = ax.hist(
        a, bins, density=1, edgecolor="skyblue", color="skyblue", alpha=0.5
    )
    mu = np.mean(a)
    sigma = np.std(a, axis=0, ddof=1).item()
    risk = [
        mu,
        mu - sigma,
        mu - rk.MAD(a),
        -rk.VaR_Hist(a, alpha),
        -rk.CVaR_Hist(a, alpha),
        -rk.EVaR_Hist(a, alpha)[0],
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
        "{0:.2%}".format((1 - alpha)) + " Confidence VaR: " + "{0:.2%}".format(risk[3]),
        "{0:.2%}".format((1 - alpha))
        + " Confidence CVaR: "
        + "{0:.2%}".format(risk[4]),
        "{0:.2%}".format((1 - alpha))
        + " Confidence EVaR: "
        + "{0:.2%}".format(risk[5]),
        "Worst Realization: " + "{0:.2%}".format(risk[6]),
    ]
    color = ["b", "r", "fuchsia", "darkorange", "limegreen", "dodgerblue", "darkgrey"]

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
        label="Normal: $\mu="
        + "{0:.2%}".format(mu)
        + "$%, $\sigma="
        + "{0:.2%}".format(sigma)
        + "$%",
    )

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

    fig = plt.gcf()
    fig.tight_layout()

    return ax


def plot_drawdown(nav, w, alpha=0.05, height=8, width=10, ax=None):
    r"""
    Create a chart with the evolution of portfolio prices and drawdown.

    Parameters
    ----------
    nav : DataFrame
        Cumulative assets returns.
    w : DataFrame, optional
        A portfolio specified by the user to compare with the efficient
        frontier. The default is None.
    alpha : float, optional
        Significante level of DaR and CDaR. The default is 0.05.
    height : float, optional
        Height of the image in inches. The default is 8.
    width : float, optional
        Width of the image in inches. The default is 10.
    ax : matplotlib axis of size (2,1), optional
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

        nav=port.nav

        ax = plf.plot_drawdown(nav=nav, w=w1, alpha=0.05, height=8, width=10, ax=None)

    .. image:: images/Drawdown.png


    """

    if not isinstance(nav, pd.DataFrame):
        raise ValueError("nav must be a DataFrame")

    if not isinstance(w, pd.DataFrame):
        raise ValueError("w must be a DataFrame")

    if w.shape[1] > 1 and w.shape[0] == 0:
        w = w.T
    elif w.shape[1] > 1 and w.shape[0] > 0:
        raise ValueError("w must be a  DataFrame")

    if nav.shape[1] != w.shape[0]:
        a1 = str(nav.shape)
        a2 = str(w.shape)
        raise ValueError("shapes " + a1 + " and " + a2 + " not aligned")

    if ax is None:
        fig = plt.gcf()
        ax = fig.subplots(nrows=2, ncols=1)
        ax = ax.flatten()
        fig.set_figwidth(width)
        fig.set_figheight(height)

    index = nav.index.tolist()

    a = np.array(nav, ndmin=2)
    a = np.insert(a, 0, 0, axis=0)
    a = np.diff(a, axis=0)
    a = np.array(a, ndmin=2) @ np.array(w, ndmin=2)
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
        "Maximum Drawdown: " + "{0:.2%}".format(risk[5]),
    ]
    color2 = ["b", "r", "fuchsia", "limegreen", "dodgerblue", "darkgrey"]

    j = 0

    ymin = np.min(DD) * 1.5

    for i in ax:
        i.clear()
        i.plot_date(index, data[j], "-", color=color1[j])
        if j == 1:
            i.fill_between(index, 0, data[j], facecolor=color1[j], alpha=0.3)
            for k in range(0, len(risk)):
                i.axhline(y=risk[k], color=color2[k], linestyle="-", label=label[k])
            i.set_ylim(ymin, 0)
            i.legend(loc="lower right")  # , fontsize = 'x-small')
        i.set_title(titles[j])
        i.xaxis.set_major_locator(
            mdates.AutoDateLocator(tz=None, minticks=5, maxticks=10)
        )
        i.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))
        ticks_loc = i.get_yticks().tolist()
        i.set_yticks(i.get_yticks())
        i.set_yticklabels(["{:3.2%}".format(x) for x in ticks_loc])
        i.grid(linestyle=":")
        j = j + 1

    fig = plt.gcf()
    fig.tight_layout()

    return ax


def plot_table(
    returns, w, MAR=0, alpha=0.05, height=9, width=12, t_factor=252, ax=None
):
    r"""
    Create a table with information about risk measures and risk adjusted
    return ratios.

    Parameters
    ----------
    returns : DataFrame
        Assets returns.
    w : DataFrame
        Portfolio weights.
    MAR: float, optional
        Minimum acceptable return.
    alpha: float, optional
        Significance level for VaR, CVaR, EVaR, DaR and CDaR.
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

        ax = plf.plot_table(returns=Y, w=w1, MAR=0, alpha=0.05, ax=None)

    .. image:: images/Port_Table.png


    """
    if not isinstance(returns, pd.DataFrame):
        raise ValueError("returns must be a DataFrame")

    if not isinstance(w, pd.DataFrame):
        raise ValueError("w must be a DataFrame")

    if returns.shape[1] != w.shape[0]:
        a1 = str(returns.shape)
        a2 = str(w.shape)
        raise ValueError("shapes " + a1 + " and " + a2 + " not aligned")

    if ax is None:
        ax = plt.gca()
        fig = plt.gcf()
        fig.set_figwidth(width)
        fig.set_figheight(height)

    mu = returns.mean()
    cov = returns.cov()
    days = (returns.index[-1] - returns.index[0]).days + 1

    X = returns @ w
    X = X.to_numpy().ravel()

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
        "Max Drawdown (MDD)",
        "(1) Annualized, multiplied by " + str(t_factor),
        "(2) Annualized, multiplied by √" + str(t_factor),
        "(3) Based on uncompounded cumulated returns",
    ]

    indicators = [
        "",
        (mu @ w).to_numpy().item() * t_factor,
        np.power(np.prod(1 + X), 360 / days) - 1,
        MAR,
        alpha,
        "",
        "",
        np.sqrt(w.T @ cov @ w).to_numpy().item() * t_factor ** 0.5,
        rk.MAD(X) * t_factor ** 0.5,
        rk.SemiDeviation(X) * t_factor ** 0.5,
        rk.LPM(X, MAR=MAR, p=1) * t_factor ** 0.5,
        rk.LPM(X, MAR=MAR, p=2) * t_factor ** 0.5,
        rk.VaR_Hist(X, alpha=alpha) * t_factor ** 0.5,
        rk.CVaR_Hist(X, alpha=alpha) * t_factor ** 0.5,
        rk.EVaR_Hist(X, alpha=alpha)[0] * t_factor ** 0.5,
        rk.WR(X) * t_factor ** 0.5,
        st.skew(X, bias=False),
        st.kurtosis(X, bias=False),
        "",
        "",
        rk.UCI_Abs(X),
        rk.ADD_Abs(X),
        rk.DaR_Abs(X),
        rk.CDaR_Abs(X, alpha=alpha),
        rk.EDaR_Abs(X, alpha=alpha)[0],
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
            if j in [6, 19]:
                cellDict[(j, 0)].set_facecolor("white")
                cellDict[(j, i)].set_facecolor("white")
            if j in [1, 7, 20]:
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

    fig = plt.gcf()
    fig.tight_layout()

    return ax
