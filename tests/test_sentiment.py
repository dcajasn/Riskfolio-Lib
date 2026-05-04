# Copyright (C) 2020-2026 Dany Cajas

import json
from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest

import riskfolio as rp


class _Response:
    def __init__(self, payload):
        self._payload = payload

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, traceback):
        return False

    def read(self):
        return json.dumps(self._payload).encode("utf-8")


def test_adanos_sentiment_downloads_compare_data():
    payload = {
        "stocks": [
            {"ticker": "AAPL", "sentiment_score": 0.35, "buzz_score": 71.2},
            {"ticker": "MSFT", "sentiment_score": -0.2, "buzz_score": 42.5},
        ]
    }

    with patch("urllib.request.urlopen", return_value=_Response(payload)) as urlopen:
        data = rp.adanos_sentiment(
            "aapl, msft",
            source="News",
            days=3,
            api_key="test-key",
            base_url="https://api.adanos.org/",
            timeout=5,
        )

    request = urlopen.call_args.args[0]
    assert request.full_url == (
        "https://api.adanos.org/news/stocks/v1/compare?tickers=AAPL%2CMSFT&days=3"
    )
    assert request.headers["X-api-key"] == "test-key"
    assert request.headers["Accept"] == "application/json"
    assert urlopen.call_args.kwargs["timeout"] == 5
    assert list(data.index) == ["AAPL", "MSFT"]
    assert data.loc["AAPL", "sentiment_score"] == 0.35


def test_adanos_sentiment_uses_environment_api_key(monkeypatch):
    monkeypatch.setenv("ADANOS_API_KEY", "env-key")

    with patch("urllib.request.urlopen", return_value=_Response({"stocks": []})) as urlopen:
        data = rp.adanos_sentiment(["AAPL"], api_key="")

    request = urlopen.call_args.args[0]
    assert request.headers["X-api-key"] == "env-key"
    assert data.empty
    assert data.index.name == "ticker"


def test_adanos_sentiment_validates_inputs(monkeypatch):
    monkeypatch.delenv("ADANOS_API_KEY", raising=False)

    with pytest.raises(ValueError, match="tickers must contain"):
        rp.adanos_sentiment([])

    with pytest.raises(ValueError, match="tickers must contain"):
        rp.adanos_sentiment(None)

    with pytest.raises(ValueError, match="source must be one of"):
        rp.adanos_sentiment(["AAPL"], source="invalid", api_key="test-key")

    with pytest.raises(ValueError, match="Please provide api_key"):
        rp.adanos_sentiment(["AAPL"], api_key="")


def test_adanos_sentiment_views_builds_black_litterman_views():
    sentiment = pd.DataFrame(
        {
            "sentiment_score": [0.5, -0.25, 0.0],
            "buzz_score": [70.0, 55.0, 10.0],
        },
        index=["AAPL", "MSFT", "TSLA"],
    )

    P, Q = rp.adanos_sentiment_views(
        sentiment,
        ["MSFT", "AAPL", "NVDA"],
        scale=0.02,
        threshold=0.1,
    )

    np.testing.assert_array_equal(P, np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]))
    np.testing.assert_array_almost_equal(Q, np.array([[-0.005], [0.01]]))


def test_adanos_sentiment_views_validates_data():
    with pytest.raises(ValueError, match="sentiment must be a DataFrame"):
        rp.adanos_sentiment_views([], ["AAPL"])

    with pytest.raises(ValueError, match="score column"):
        rp.adanos_sentiment_views(pd.DataFrame(index=["AAPL"]), ["AAPL"])

    sentiment = pd.DataFrame({"sentiment_score": [0.2, 0.3]}, index=["AAPL", "AAPL"])
    with pytest.raises(ValueError, match="duplicate tickers"):
        rp.adanos_sentiment_views(sentiment, ["AAPL"])


def test_adanos_sentiment_views_returns_empty_view_matrices():
    sentiment = pd.DataFrame({"sentiment_score": [0.05]}, index=["AAPL"])

    P, Q = rp.adanos_sentiment_views(sentiment, ["AAPL", "MSFT"], threshold=0.1)

    assert P.shape == (0, 2)
    assert Q.shape == (0, 1)
