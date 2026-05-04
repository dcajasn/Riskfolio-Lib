""""""  #

"""
Copyright (c) 2020-2026, Dany Cajas
All rights reserved.
This work is licensed under BSD 3-Clause "New" or "Revised" License.
License available at https://github.com/dcajasn/Riskfolio-Lib/blob/master/LICENSE.txt
"""

import json
import os
import urllib.parse
import urllib.request

import numpy as np
import pandas as pd

__all__ = [
    "adanos_sentiment",
    "adanos_sentiment_views",
]

ADANOS_SOURCES = {"reddit", "x", "news", "polymarket"}
ADANOS_BASE_URL = "https://api.adanos.org"


def adanos_sentiment(
    tickers,
    source="reddit",
    days=7,
    api_key=None,
    base_url=ADANOS_BASE_URL,
    timeout=10,
):
    r"""
    Download Adanos Market Sentiment API data for US stocks.

    Parameters
    ----------
    tickers : str or iterable
        Tickers to query. It can be a comma separated string like
        ``"AAPL,MSFT"`` or an iterable like ``["AAPL", "MSFT"]``.

    source : str, optional
        Adanos source. Possible values are ``"reddit"``, ``"x"``,
        ``"news"`` and ``"polymarket"``. The default is ``"reddit"``.

    days : int, optional
        Lookback window in days. The default is 7.

    api_key : str, optional
        Adanos API key. If not provided, ``ADANOS_API_KEY`` environment
        variable will be used.

    base_url : str, optional
        Adanos API base URL. The default is ``"https://api.adanos.org"``.

    timeout : int or float, optional
        Request timeout in seconds. The default is 10.

    Returns
    -------
    data : DataFrame
        Sentiment data indexed by ticker.

    Examples
    --------
    ::

        import riskfolio as rp

        sentiment = rp.adanos_sentiment(["AAPL", "MSFT"], source="news")

    """

    source = _normalize_source(source)
    ticker_list = _normalize_tickers(tickers)
    api_key = _get_api_key(api_key)
    params = urllib.parse.urlencode({"tickers": ",".join(ticker_list), "days": days})
    url = "{}/{}/stocks/v1/compare?{}".format(base_url.rstrip("/"), source, params)
    request = urllib.request.Request(
        url,
        headers={"X-API-Key": api_key, "Accept": "application/json"},
        method="GET",
    )

    with urllib.request.urlopen(request, timeout=timeout) as response:
        payload = json.loads(response.read().decode("utf-8"))

    data = pd.DataFrame(payload.get("stocks", []))
    if data.empty:
        return pd.DataFrame(index=pd.Index([], name="ticker"))
    if "ticker" not in data.columns:
        raise ValueError("Adanos response must include a ticker field")

    data["ticker"] = data["ticker"].astype(str).str.upper()
    data = data.set_index("ticker")
    data.index.name = "ticker"
    return data


def adanos_sentiment_views(
    sentiment,
    assets,
    score="sentiment_score",
    scale=0.01,
    threshold=0.0,
):
    r"""
    Build Black Litterman absolute views from Adanos sentiment scores.

    Parameters
    ----------
    sentiment : DataFrame
        Sentiment data indexed by ticker. It can be the output of
        :func:`adanos_sentiment`.

    assets : list
        Assets in the same order used by the portfolio returns matrix.

    score : str, optional
        Column used to build the return views. The default is
        ``"sentiment_score"``.

    scale : scalar, optional
        Multiplier applied to the score to convert sentiment units into
        return views. The default is 0.01.

    threshold : scalar, optional
        Minimum absolute score required to include a view. The default is 0.

    Returns
    -------
    P : nd-array
        Matrix of shape (n_views, n_assets) that identifies the assets in
        each view.

    Q : nd-array
        Expected return vector of shape (n_views, 1) implied by the sentiment
        views.

    Examples
    --------
    ::

        import riskfolio as rp

        sentiment = rp.adanos_sentiment(["AAPL", "MSFT"], source="news")
        P, Q = rp.adanos_sentiment_views(sentiment, ["AAPL", "MSFT"], scale=0.02)

    """

    if not isinstance(sentiment, pd.DataFrame):
        raise ValueError("sentiment must be a DataFrame")
    if score not in sentiment.columns:
        raise ValueError("sentiment must include the score column")

    assets = _normalize_tickers(assets)
    sentiment = sentiment.copy()
    sentiment.index = sentiment.index.astype(str).str.upper()
    if sentiment.index.has_duplicates:
        raise ValueError("sentiment index must not contain duplicate tickers")

    P = []
    Q = []
    for i, asset in enumerate(assets):
        if asset not in sentiment.index:
            continue
        value = sentiment.loc[asset, score]
        if pd.isna(value):
            continue
        value = float(value)
        if abs(value) < threshold:
            continue
        view = [0.0] * len(assets)
        view[i] = 1.0
        P.append(view)
        Q.append([value * scale])

    if len(P) == 0:
        return np.empty((0, len(assets))), np.empty((0, 1))

    return np.array(P, ndmin=2, dtype=float), np.array(Q, ndmin=2, dtype=float)


def _get_api_key(api_key):
    api_key = (api_key or os.getenv("ADANOS_API_KEY", "")).strip()
    if not api_key:
        raise ValueError("Please provide api_key or set ADANOS_API_KEY")
    return api_key


def _normalize_source(source):
    source = source.strip().lower()
    if source not in ADANOS_SOURCES:
        raise ValueError("source must be one of: reddit, x, news, polymarket")
    return source


def _normalize_tickers(tickers):
    if tickers is None:
        raise ValueError("tickers must contain at least one ticker")
    if isinstance(tickers, str):
        tickers = tickers.split(",")
    tickers = [str(ticker).strip().upper() for ticker in tickers if str(ticker).strip()]
    if not tickers:
        raise ValueError("tickers must contain at least one ticker")
    return tickers
