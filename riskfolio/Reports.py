import pandas as pd
import matplotlib.pyplot as plt
from xlsxwriter.utility import xl_range_abs, xl_range, xl_rowcol_to_cell, xl_col_to_name
import datetime
import riskfolio.PlotFunctions as plf
import riskfolio.RiskFunctions as rk


__LICENSE__ = """Copyright (c) 2020-2021, Dany Cajas
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
* Redistributions of source code must retain the above copyright notice, this  list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
* Neither the name of Riskfolio-Lib nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."""


def jupyter_report(
    returns,
    w,
    rm="MV",
    rf=0,
    alpha=0.05,
    others=0.05,
    nrow=25,
    height=6,
    width=14,
    t_factor=252,
):
    r"""
    Create a matplotlib report with useful information to analyze risk and
    profitability of investment portfolios.

    Parameters
    ----------
    returns : DataFrame
        Assets returns.
    w : DataFrame of shape (n_assets, 1)
        Portfolio weights.
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
        - 'EDaR': Entropic Drawdown at Risk of uncompounded cumulative returns.
        - 'UCI': Ulcer Index of uncompounded cumulative returns.

    rf : float, optional
        Risk free rate or minimum aceptable return. The default is 0.
    alpha : float, optional
        Significante level of VaR, CVaR, EVaR, DaR and CDaR.
        The default is 0.05.
    others : float, optional
        Percentage of others section. The default is 0.05.
    nrow : int, optional
        Number of rows of the legend. The default is 25.
    height : float, optional
        Average height of charts in the image in inches. The default is 6.
    width : float, optional
        Width of the image in inches. The default is 14.
    t_factor : float, optional
        Factor used to annualize expected return and expected risks for
        risk measures based on returns (not drawdowns). The default is 252.
        
        .. math::
            
            \begin{align}
            \text{Annualized Return} & = \text{Return} \, \times \, \text{t_factor} \\
            \text{Annualized Risk} & = \text{Risk} \, \times \, \sqrt{\text{t_factor}}
            \end{align}
        
    ax : matplotlib axis of size (6,1), optional
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

        ax = rp.jupyter_report(returns, w, rm='MV', rf=0, alpha=0.05, height=6, width=14,
                               others=0.05, nrow=25)

    .. image:: images/Report_1.png
    .. image:: images/Report_2.png
    .. image:: images/Report_3.png
    .. image:: images/Report_4.png

    """

    cov = returns.cov()
    nav = returns.cumsum()

    fig, ax = plt.subplots(
        nrows=6,
        figsize=(width, height * 6),
        gridspec_kw={"height_ratios": [2, 1, 1.5, 1, 1, 1]},
    )

    ax[0] = plf.plot_table(returns, w, MAR=rf, alpha=alpha, t_factor=t_factor, ax=ax[0])

    ax[2] = plf.plot_pie(
        w=w,
        title="Portfolio Composition",
        others=others,
        nrow=nrow,
        cmap="tab20",
        ax=ax[2],
    )

    ax[3] = plf.plot_risk_con(
        w=w, cov=cov, returns=returns, rm=rm, rf=rf, alpha=alpha, ax=ax[3]
    )

    ax[4] = plf.plot_hist(returns=returns, w=w, alpha=alpha, bins=50, ax=ax[4])

    ax[[1, 5]] = plf.plot_drawdown(nav=nav, w=w, alpha=0.05, ax=ax[[1, 5]])

    year = str(datetime.datetime.now().year)

    title = "Riskfolio-Lib Report"
    subtitle = "Copyright (c) 2020-" + year + ", Dany Cajas. All rights reserved."

    fig.suptitle(title, fontsize="xx-large", y=1.011, fontweight="bold")
    ax[0].set_title(subtitle, fontsize="large", ha="center", pad=10)

    return ax


def excel_report(returns, w, rf=0, alpha=0.05, t_factor=252, name="report"):
    r"""
    Create an Excel report (with formulas) with useful information to analyze
    risk and profitability of investment portfolios.

    Parameters
    ----------
    returns : DataFrame
        Assets returns.
    w : DataFrame of size (n_assets, n_portfolios)
        Portfolio weights.
    rf : float, optional
        Risk free rate or minimum aceptable return. The default is 0.
    alpha : float, optional
        Significante level of VaR, CVaR, EVaR, DaR and CDaR.
        The default is 0.05.
    t_factor : float, optional
        Factor used to annualize expected return and expected risks for
        risk measures based on returns (not drawdowns). The default is 252.
        
        .. math::
            
            \begin{align}
            \text{Annualized Return} & = \text{Return} \, \times \, \text{t_factor} \\
            \text{Annualized Risk} & = \text{Risk} \, \times \, \sqrt{\text{t_factor}}
            \end{align}
        
    name : str, optional
        Name or name with path where the Excel report will be saved. If no
        path is provided the report will be saved in the same path of
        current file.

    Raises
    ------
    ValueError
        When the report cannot be built.

    Example
    -------
    ::

        rp.excel_report(returns, w, MAR=0, alpha=0.05, name='report', files=None)

    .. image:: images/Excel.png

    """
    n1 = w.shape[0]
    n2 = returns.shape[0]

    portfolios = w.columns.tolist()
    dates = returns.index.tolist()
    year = str(datetime.datetime.now().year)
    days = (returns.index[-1] - returns.index[0]).days + 1

    # Create a Pandas Excel writer using XlsxWriter as the engine.
    writer = pd.ExcelWriter(name + ".xlsx", engine="xlsxwriter")

    # Convert the dataframe to an XlsxWriter Excel object.
    w.to_excel(writer, sheet_name="Resume", startrow=36, startcol=0)
    returns.to_excel(writer, sheet_name="Returns", index_label=["Date"])

    # Get the xlsxwriter objects from the dataframe writer object.
    workbook = writer.book
    worksheet1 = writer.sheets["Resume"]
    worksheet2 = writer.sheets["Returns"]
    worksheet3 = workbook.add_worksheet("Portfolios")
    worksheet4 = workbook.add_worksheet("Absdev")
    worksheet5 = workbook.add_worksheet("CumRet")
    worksheet6 = workbook.add_worksheet("Drawdown")
    worksheet7 = workbook.add_worksheet("devBelowTarget")
    worksheet8 = workbook.add_worksheet("devBelowMean")

    worksheet1.hide_gridlines(2)
    worksheet2.hide_gridlines(2)
    worksheet3.hide_gridlines(2)
    worksheet4.hide_gridlines(2)
    worksheet5.hide_gridlines(2)
    worksheet6.hide_gridlines(2)
    worksheet7.hide_gridlines(2)
    worksheet8.hide_gridlines(2)

    # Cell Formats
    cell_format1 = workbook.add_format({"bold": True, "border": True})
    cell_format2 = workbook.add_format({"bold": True, "font_size": 28, "right": True})
    cell_format3 = workbook.add_format({"num_format": "0.0000%"})
    cell_format4 = workbook.add_format({"num_format": "0.0000%", "border": True})
    cell_format5 = workbook.add_format({"num_format": "yyyy-mm-dd", "bold": True})
    cell_format6 = workbook.add_format({"num_format": "0.0000", "border": True})
    cell_format7 = workbook.add_format(
        {"num_format": "yyyy-mm-dd", "bold": True, "border": True}
    )
    cell_format8 = workbook.add_format({"num_format": "0,000", "border": True})

    cols = xl_col_to_name(1) + ":" + xl_col_to_name(n2)
    worksheet1.set_column(cols, 11, cell_format3)
    worksheet2.set_column(cols, 9, cell_format3)

    worksheet2.write(0, 0, "Date", cell_format1)
    worksheet3.write(0, 0, "Date", cell_format1)
    worksheet4.write(0, 0, "Date", cell_format1)
    worksheet5.write(0, 0, "Date", cell_format1)
    worksheet6.write(0, 0, "Date", cell_format1)
    worksheet7.write(0, 0, "Date", cell_format1)
    worksheet8.write(0, 0, "Date", cell_format1)

    worksheet1.set_column("A:A", 36)
    worksheet2.set_column("A:A", 10, cell_format5)
    worksheet3.set_column("A:A", 10, cell_format5)
    worksheet4.set_column("A:A", 10, cell_format5)
    worksheet5.set_column("A:A", 10, cell_format5)
    worksheet6.set_column("A:A", 10, cell_format5)
    worksheet7.set_column("A:A", 10, cell_format5)
    worksheet8.set_column("A:A", 10, cell_format5)

    for i in range(0, n2):
        r = xl_rowcol_to_cell(i + 1, 0)
        formula = "=Returns!" + r + ""
        worksheet2.write(i + 1, 0, dates[i], cell_format7)
        worksheet3.write_formula(i + 1, 0, formula, cell_format7)
        worksheet4.write_formula(i + 1, 0, formula, cell_format7)
        worksheet5.write_formula(i + 1, 0, formula, cell_format7)
        worksheet6.write_formula(i + 1, 0, formula, cell_format7)
        worksheet7.write_formula(i + 1, 0, formula, cell_format7)
        worksheet8.write_formula(i + 1, 0, formula, cell_format7)

    labels_1 = [
        "",
        "",
        "",
        "",
        "Profitability and Other Inputs",
        "Total Days in DataBase",
        "Mean Return (1)",
        "Compound Annual Growth Rate (CAGR)",
        "Minimum Acceptable Return (MAR) (1)",
        "Alpha",
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
        "Entropic Drawdown at Risk (CDaR)",
        "Max Drawdown (MDD)",
    ]

    for i in range(0, len(labels_1)):
        if labels_1[i] != "":
            worksheet1.write(i, 0, labels_1[i], cell_format1)

    for i in range(0, len(portfolios)):
        a = "Portfolio " + str(i + 1)
        worksheet1.write(3, 1 + i, a, cell_format1)
        worksheet1.write(36, 1 + i, a, cell_format1)
        worksheet3.write(0, 1 + i, a, cell_format1)
        worksheet4.write(0, 1 + i, a, cell_format1)
        worksheet5.write(0, 1 + i, a, cell_format1)
        worksheet6.write(0, 1 + i, a, cell_format1)
        worksheet7.write(0, 1 + i, a, cell_format1)
        worksheet8.write(0, 1 + i, a, cell_format1)

    for j in range(0, len(portfolios)):
        r_0 = xl_rowcol_to_cell(8, 1 + j)  # MAR cell
        r_1 = xl_range_abs(37, 1 + j, 36 + n1, 1 + j)
        r_2 = xl_range_abs(1, 1 + j, n2, 1 + j)
        for i in range(0, n2):
            r_3 = xl_range(i + 1, 1, i + 1, n1)
            r_4 = xl_rowcol_to_cell(i + 1, 1 + j)
            r_5 = xl_range_abs(1, 1 + j, i + 1, 1 + j)
            formula1 = "{=MMULT(" + "Returns!" + r_3 + ",Resume!" + r_1 + ")}"
            formula2 = "=ABS(Portfolios!" + r_4 + "-AVERAGE(Portfolios!" + r_2 + "))"
            formula3 = "=SUM(Portfolios!" + r_5 + ")"
            formula4 = "=MAX(CumRet!" + r_5 + ")-CumRet!" + r_4
            formula5 = (
                "=MAX(Resume!"
                + r_0
                + "/ "
                + str(t_factor)
                + "-Portfolios!"
                + r_4
                + ", 0)"
            )
            formula6 = "=MAX(AVERAGE(Portfolios!" + r_2 + ")-Portfolios!" + r_4 + ", 0)"
            worksheet3.write_formula(i + 1, 1 + j, formula1, cell_format3)
            worksheet4.write_formula(i + 1, 1 + j, formula2, cell_format3)
            worksheet5.write_formula(i + 1, 1 + j, formula3, cell_format3)
            worksheet6.write_formula(i + 1, 1 + j, formula4, cell_format3)
            worksheet7.write_formula(i + 1, 1 + j, formula5, cell_format3)
            worksheet8.write_formula(i + 1, 1 + j, formula6, cell_format3)

        r_6 = xl_rowcol_to_cell(9, 1 + j)  # Alpha cell
        r_7 = xl_rowcol_to_cell(17, 1 + j)  # Value at Risk cell
        AVG = "=AVERAGE(Portfolios!" + r_2 + ") * " + str(t_factor) + ""
        CUM = "{=PRODUCT(1 + Portfolios!" + r_2 + ")^(360/" + str(days) + ")-1}"
        STDEV = "=STDEV(Portfolios!" + r_2 + ") * SQRT(" + str(t_factor) + ")"
        MAD = "=AVERAGE(Absdev!" + r_2 + ") * SQRT(" + str(t_factor) + ")"
        ALPHA = "=" + str(alpha)
        VaR = (
            "=-SMALL(Portfolios!"
            + r_2
            + ",ROUNDUP(COUNT(Portfolios!"
            + r_2
            + ")*"
            + r_6
            + ",0)) * SQRT("
            + str(t_factor)
            + ")"
        )
        CVaR = (
            "=-((SUMIF(Portfolios!"
            + r_2
            + ',"<="&(-'
            + r_7
            + "/SQRT("
            + str(t_factor)
            + ")),Portfolios!"
            + r_2
            + ")"
        )
        CVaR += (
            "-ROUNDUP(COUNT(Portfolios!"
            + r_2
            + ")*"
            + r_6
            + ",0)*(-"
            + r_7
            + "/SQRT("
            + str(t_factor)
            + ")))/(COUNT(Portfolios!"
            + r_2
            + ")*"
            + r_6
            + ")-"
            + r_7
            + "/SQRT("
            + str(t_factor)
            + ")) * SQRT("
            + str(t_factor)
            + ")"
        )
        EVaR = (
            "="
            + str(rk.EVaR_Hist(returns @ w, alpha=alpha)[0])
            + " * SQRT("
            + str(t_factor)
            + ")"
        )
        WR = "=-MIN(Portfolios!" + r_2 + ") * SQRT(" + str(t_factor) + ")"
        MDD = "=MAX(Drawdown!" + r_2 + ")"
        ADD = "=AVERAGE(Drawdown!" + r_2 + ")"
        DaR = (
            "=+LARGE(Drawdown!"
            + r_2
            + ",ROUNDUP(COUNT(Drawdown!"
            + r_2
            + ")*"
            + r_6
            + ",0))"
        )
        CDaR = (
            "=((SUMIF(Drawdown!" + r_2 + ',">="&' + DaR[2:] + ",Drawdown!" + r_2 + ")"
        )
        CDaR += (
            "-ROUNDUP(COUNT(Drawdown!"
            + r_2
            + ")*"
            + r_6
            + ",0)*"
            + DaR[2:]
            + ")/(COUNT(Drawdown!"
            + r_2
            + ")*"
            + r_6
            + ")+"
            + DaR[2:]
            + ")"
        )
        EDaR = "=" + str(rk.EDaR_Abs(returns @ w, alpha=alpha)[0])
        UCI = "=SQRT(SUMSQ(Drawdown!" + r_2 + ")/COUNT(Drawdown!" + r_2 + "))"
        MAR = "=" + str(rf)
        FLPM = "=AVERAGE(devBelowTarget!" + r_2 + ") * SQRT(" + str(t_factor) + ")"
        SLPM = (
            "=SQRT(SUMSQ(devBelowTarget!"
            + r_2
            + ")/(COUNT(devBelowTarget!"
            + r_2
            + ") - 1))"
            + " * SQRT("
            + str(t_factor)
            + ")"
        )
        SDEV = (
            "=SQRT(SUMSQ(devBelowMean!"
            + r_2
            + ")/(COUNT(devBelowMean!"
            + r_2
            + ") - 1))"
            + " * SQRT("
            + str(t_factor)
            + ")"
        )
        SKEW = "=SKEW(Portfolios!" + r_2 + ")"
        KURT = "=KURT(Portfolios!" + r_2 + ")"

        labels_2 = [
            "",
            "",
            "",
            "",
            "",
            str(days),
            AVG,
            CUM,
            MAR,
            ALPHA,
            "",
            "",
            STDEV,
            MAD,
            SDEV,
            FLPM,
            SLPM,
            VaR,
            CVaR,
            EVaR,
            WR,
            SKEW,
            KURT,
            "",
            "",
            UCI,
            ADD,
            DaR,
            CDaR,
            EDaR,
            MDD,
        ]

        for i in range(0, len(labels_2)):
            if labels_1[i] in ["Skewness", "Kurtosis"]:
                worksheet1.write_formula(i, 1 + j, labels_2[i], cell_format6)
            elif labels_1[i] in ["Total Days in DataBase"]:
                worksheet1.write_formula(i, 1 + j, labels_2[i], cell_format8)
            elif labels_2[i] != "":
                worksheet1.write_formula(i, 1 + j, labels_2[i], cell_format4)

    merge_format = workbook.add_format({"align": "Left", "valign": "vjustify"})
    merge_format.set_text_wrap()
    worksheet1.set_row(1, 215)
    worksheet1.merge_range("A2:K2", __LICENSE__.replace("2021", year), merge_format)
    worksheet1.write(31, 0, "(1) Annualized, multiplied by " + str(t_factor))
    worksheet1.write(32, 0, "(2) Annualized, multiplied by √" + str(t_factor))
    worksheet1.write(33, 0, "(3) Based on uncompounded cumulated returns")
    worksheet1.write(0, 0, "Riskfolio-Lib Report", cell_format2)

    writer.save()
    workbook.close()
