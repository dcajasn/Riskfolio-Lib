import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
import riskfolio.RiskFunctions as rk

__all__ = [
    "plot_series",
    "plot_frontier",
    "plot_pie",
    "plot_frontier_area",
    "plot_hist",
    "plot_drawdown",
]

rm_names = [
    "Standard Deviation",
    "Mean Absolute Deviation",
    "Semi Standard Deviation",
    "Value at Risk",
    "Conditional Value at Risk",
    "Worst Realization",
    "First Lower Partial Moment",
    "Second Lower Partial Moment",
    "Max Drawdown",
    "Average Drawdown",
    "Conditional Drawdown at Risk",
]

rmeasures = [
    "MV",
    "MAD",
    "MSV",
    "VaR",
    "CVaR",
    "WR",
    "FLPM",
    "SLPM",
    "MDD",
    "ADD",
    "CDaR",
]


def plot_series(returns, w, cmap="tab20", height=6, width=10, ax=None):
    r"""
    Create a chart with the compound cumulated of the portfolios.

    Parameters
    ----------
    returns : DataFrame
        Assets returns.
    w : DataFrame
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
        
        ax = plf.plot_series(data=Y, w=ws, cmap='tab20', height=6, width=10, ax=None)
        
    .. image:: images/Port_Series.png


    """
    if not isinstance(returns, pd.DataFrame):
        raise ValueError("data must be a DataFrame")

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
        a = np.matrix(returns) * np.matrix(w[X[i]]).T
        prices = 1 + np.insert(a, 0, 0, axis=0)
        prices = np.cumprod(prices, axis=0)
        prices = np.ravel(prices).tolist()
        del prices[0]

        ax.plot_date(index, prices, "-", label=labels[i])

    ax.set_yticklabels(["{:3.2f}".format(x) for x in ax.get_yticks()])
    ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))

    fig.tight_layout()

    return ax


def plot_frontier(
    w_frontier,
    mu,
    cov=None,
    returns=None,
    rm="MV",
    rf=0,
    alpha=0.01,
    cmap="viridis",
    w=None,
    label="Portfolio",
    marker="*",
    s=16,
    c="r",
    height=6,
    width=10,
    ax=None,
):
    """
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
        Risk measure used to create the frontier. The default is 'MV'.
    rf : float, optional
        Risk free rate or minimum aceptable return. The default is 0.
    alpha : float, optional
        Significante level of VaR, CVaR and CDaR. The default is 0.01.
    cmap : cmap, optional
        Colorscale, represente the risk adjusted return ratio.
        The default is 'viridis'.
    w : DataFrame, optional
        A portfolio specified by the user. The default is None.
    label : str, optional
        Name of portfolio that appear on plot legend.
        The default is 'Portfolio'.
    marker : str, optional
        Marker of w_. The default is '*'.
    s : float, optional
        Size of marker. The default is 16.
    c : str, optional
        Color of marker. The default is 'r'.
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
                               rm=rm, rf=0, alpha=0.01, cmap='viridis', w=w1,
                               label='Portfolio', marker='*', s=16, c='r',
                               height=6, width=10, ax=None)
        
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

    mu_ = np.matrix(mu)

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
        weights = np.matrix(w_frontier.iloc[:, i]).T
        risk = rk.Sharpe_Risk(
            weights, cov=cov, returns=returns, rm=rm, rf=rf, alpha=alpha
        )
        ret = mu_ * weights
        ret = ret.item()
        ratio = (ret - rf) / risk

        X1.append(risk)
        Y1.append(ret)
        Z1.append(ratio)

    ax1 = ax.scatter(X1, Y1, c=Z1, cmap=cmap)

    if w is not None:
        X2 = []
        Y2 = []
        for i in range(w.shape[1]):
            weights = np.matrix(w.iloc[:, i]).T
            risk = rk.Sharpe_Risk(
                weights, cov=cov, returns=returns, rm=rm, rf=rf, alpha=alpha
            )
            ret = mu_ * weights
            ret = ret.item()
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

    ax.set_yticklabels(["{:.4%}".format(x) for x in ax.get_yticks()])
    ax.set_xticklabels(["{:.4%}".format(x) for x in ax.get_xticks()])

    ax.tick_params(axis="y", direction="in")
    ax.tick_params(axis="x", direction="in")

    ax.grid(linestyle=":")

    colorbar = ax.figure.colorbar(ax1)
    colorbar.set_label("Risk Adjusted Return Ratio")

    fig.tight_layout()

    return ax


def plot_pie(
    w, title="", others=0.05, nrow=25, cmap="tab20", height=6, width=8, ax=None
):
    """
    Create a pie chart with portfolio weights.
    
    Parameters
    ----------
    w : DataFrame
        Weights of the portfolio.
    title : str, optional
        Title of the chart. The default is ''.
    others : float, optional
        Percentage of others section. The default is 0.05.
    nrow : int, optional
        Number of rows of the legend. The default is 25.
    cmap : cmap, optional
        Color scale, represente the risk adjusted return ratio.
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
        
        ax = plf.plot_pie(w=w1, title='Portafolio', height=6, width=10, cmap="tab20", ax=None)
        
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

    if title == "":
        title = "Portfolio Composition"

    ax.set_title(title)

    labels = w.index.tolist()
    sizes = w.iloc[:, 0].tolist()
    sizes2 = pd.DataFrame([labels, sizes]).T
    sizes2.columns = ["labels", "values"]
    sizes2 = sizes2.sort_values(by=["values"], ascending=False)
    sizes2.index = [i for i in range(0, len(labels))]
    sizes3 = sizes2.cumsum()
    l = sizes3[sizes3["values"] >= 1 - others].index.tolist()[0]

    item = pd.DataFrame(["Others", 1 - sizes3.loc[l, "values"]]).T
    item.columns = ["labels", "values"]
    sizes2 = sizes2[sizes2.index <= l]
    sizes2 = sizes2.append(item)

    sizes = sizes2["values"].tolist()
    labels = sizes2["labels"].tolist()
    sizes2 = ["{0:.1%}".format(i) for i in sizes]

    colormap = cm.get_cmap(cmap)
    colormap = colormap(np.linspace(0, 1, 20))

    if cmap == "gist_rainbow":
        colormap = colormap[::-1]

    cycle = plt.cycler("color", colormap)
    ax.set_prop_cycle(cycle)

    size = 0.4

    # set up style cycles

    wedges, texts = ax.pie(
        sizes, radius=1, wedgeprops=dict(width=size, edgecolor="black"), startangle=-15
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
        name = str(labels[i]) + " - " + str(sizes2[i])
        ax.annotate(
            name,
            xy=(x, y),
            xytext=(1.1 * np.sign(x), 1.1 * y),
            horizontalalignment=horizontalalignment,
            **kw
        )

    fig.tight_layout()

    return ax


def plot_frontier_area(w_frontier, nrow=25, cmap="tab20", height=6, width=10, ax=None):
    r"""
    Create a chart with te asset structure of the efficient frontier.
    
    Parameters
    ----------
    w_frontier : DataFrame
        Weights of portfolios in the efficient frontier.
    nrow : int, optional
        Number of rows of the legend. The default is 25.
    cmap : cmap, optional
        Color scale, represente the risk adjusted return ratio.
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
        
        ax = plf.plot_frontier_area(w_frontier=ws, cmap="tab20", height=6, width=10, ax=None)
        
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

    ax.set_yticklabels(["{:3.2%}".format(x) for x in ax.get_yticks()])
    ax.grid(linestyle=":")

    n = int(np.ceil(len(labels) / nrow))

    ax.legend(labels, loc="center left", bbox_to_anchor=(1, 0.5), ncol=n)

    fig.tight_layout()

    return ax


def plot_hist(returns, w, alpha=0.01, bins=50, height=6, width=10, ax=None):
    r"""
    Create a histogram of portfolio returns with the risk measures.
    
    Parameters
    ----------
    returns : DataFrame
        Assets returns.
    w : DataFrame, optional
        A portfolio specified by the user to compare with the efficient
        frontier. The default is None.
    alpha : float, optional
        Significante level of VaR, CVaR and CDaR. The default is 0.01.
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
        
        ax = plf.plot_hist(data=Y, w=w1, alpha=0.01, bins=50, height=6, width=10, ax=None)
        
    .. image:: images/Histogram.png
    
    """

    if not isinstance(returns, pd.DataFrame):
        raise ValueError("data must be a DataFrame")

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

    a = np.matrix(returns) * np.matrix(w)
    ax.set_title("Portfolio Returns Histogram")
    n, bins1, patches = ax.hist(
        a, bins, density=1, edgecolor="skyblue", color="skyblue", alpha=0.5
    )
    mu = np.mean(a)
    sigma = np.asscalar(np.std(a, axis=0, ddof=1))
    risk = [
        mu,
        mu - sigma,
        mu - rk.MAD(a),
        -rk.VaR_Hist(a, alpha),
        -rk.CVaR_Hist(a, alpha),
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
        "{0:.2%}".format((1 - alpha))
        + " Confidence VaR: "
        + "{0:.2%}".format(-risk[3]),
        "{0:.2%}".format((1 - alpha))
        + " Confidence CVaR: "
        + "{0:.2%}".format(-risk[4]),
        "Worst Realization: " + "{0:.2%}".format(-risk[5]),
    ]
    color = ["b", "r", "fuchsia", "darkorange", "limegreen", "darkgrey"]

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

    ax.set_xticklabels(["{:3.2%}".format(x) for x in ax.get_xticks()])
    ax.set_yticklabels(["{:3.2%}".format(x * factor) for x in ax.get_yticks()])
    ax.legend(loc="upper right")  # , fontsize = 'x-small')
    ax.grid(linestyle=":")
    ax.set_ylabel("Probability Density")

    fig.tight_layout()

    return ax


def plot_drawdown(nav, w, alpha=0.01, height=8, width=10, ax=None):
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
        Significante level of VaR, CVaR and CDaR. The default is 0.01.
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
    ax : matplotlib axis.
        Returns the Axes object with the plot for further tweaking.
    
    Example
    -------
    ::
        nav=port.nav
        
        ax = plf.plot_drawdown(nav=nav, w=w1, alpha=0.01, height=8, width=10, ax=None)

    .. image:: images/Drawdown.png
    
    """

    if not isinstance(nav, pd.DataFrame):
        raise ValueError("data must be a DataFrame")

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

    a = np.matrix(nav)
    a = np.insert(a, 0, 0, axis=0)
    a = np.diff(a, axis=0)
    a = np.matrix(a) * np.matrix(w)
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
    risk = [-rk.MaxAbsDD(a), -rk.AvgAbsDD(a), -rk.ConAbsDD(a, alpha)]
    label = [
        "Maximum Drawdown: " + "{0:.2%}".format(risk[0]),
        "Average Drawdown: " + "{0:.2%}".format(risk[1]),
        "{0:.2%}".format((1 - alpha))
        + " Confidence CDaR: "
        + "{0:.2%}".format(risk[2]),
    ]
    color2 = ["r", "limegreen", "fuchsia"]

    j = 0

    ymin = np.min(DD) * 1.4

    for i in ax:
        i.clear()
        i.plot_date(index, data[j], "-", color=color1[j])
        if j == 1:
            i.fill_between(index, 0, data[j], facecolor=color1[j], alpha=0.3)
            for k in range(0, 3):
                i.axhline(y=risk[k], color=color2[k], linestyle="-", label=label[k])
            i.set_ylim(ymin, 0)
            i.legend(loc="lower right")  # , fontsize = 'x-small')
        i.set_title(titles[j])
        i.set_yticklabels(["{:3.2%}".format(x) for x in i.get_yticks()])
        i.grid(linestyle=":")
        j = j + 1

    fig.tight_layout()

    return ax
