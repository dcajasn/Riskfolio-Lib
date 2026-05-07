# Copyright (C) 2020-2026 Dany Cajas

import ast
import json
import os
import urllib.parse
import urllib.request
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest


NOTEBOOK = (
    Path(__file__).resolve().parents[1]
    / "examples"
    / "Tutorial 54 - Black Litterman with Adanos Market Sentiment Views.ipynb"
)


class _Response:
    def __init__(self, payload):
        self._payload = payload

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, traceback):
        return False

    def read(self):
        return json.dumps(self._payload).encode("utf-8")


def _notebook():
    return json.loads(NOTEBOOK.read_text(encoding="utf-8"))


def _example_namespace():
    namespace = {
        "json": json,
        "np": np,
        "os": os,
        "pd": pd,
        "urllib": urllib,
    }
    for cell in _notebook()["cells"]:
        source = "".join(cell.get("source", []))
        if source.startswith("def adanos_compare_sentiment") or source.startswith(
            "def sentiment_to_black_litterman_views"
        ):
            source = source.split("\nP, Q = ", 1)[0]
            exec(compile(source, str(NOTEBOOK), "exec"), namespace)
    return namespace


def test_adanos_sentiment_example_notebook_code_cells_compile():
    for index, cell in enumerate(_notebook()["cells"]):
        if cell["cell_type"] == "code":
            ast.parse("".join(cell["source"]), filename=f"{NOTEBOOK}:cell-{index}")


def test_adanos_compare_sentiment_builds_expected_request():
    namespace = _example_namespace()
    payload = {
        "stocks": [
            {"ticker": "aapl", "sentiment_score": 0.25, "buzz_score": 42.0},
            {"ticker": "msft", "sentiment_score": -0.1, "buzz_score": 18.0},
        ]
    }

    with patch("urllib.request.urlopen", return_value=_Response(payload)) as urlopen:
        sentiment = namespace["adanos_compare_sentiment"](
            ["aapl", "MSFT"],
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
    assert urlopen.call_args.kwargs["timeout"] == 5
    assert list(sentiment.index) == ["AAPL", "MSFT"]
    assert sentiment.loc["AAPL", "sentiment_score"] == 0.25


def test_sentiment_to_black_litterman_views_handles_thresholds_and_shapes():
    namespace = _example_namespace()
    sentiment = pd.DataFrame(
        {"sentiment_score": [0.5, -0.25, 0.05]},
        index=["AAPL", "MSFT", "NVDA"],
    )

    P, Q = namespace["sentiment_to_black_litterman_views"](
        sentiment,
        ["MSFT", "AAPL", "NVDA"],
        annual_return_scale=0.02,
        min_abs_score=0.1,
    )

    np.testing.assert_array_equal(
        P.values,
        np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
    )
    np.testing.assert_array_almost_equal(Q.values, np.array([[-0.005], [0.01]]))
    assert list(P.index) == ["MSFT", "AAPL"]
    assert list(Q.columns) == ["views"]


def test_sentiment_to_black_litterman_views_returns_empty_view_matrices():
    namespace = _example_namespace()
    sentiment = pd.DataFrame({"sentiment_score": [0.05]}, index=["AAPL"])

    P, Q = namespace["sentiment_to_black_litterman_views"](
        sentiment,
        ["AAPL", "MSFT"],
        min_abs_score=0.1,
    )

    assert P.shape == (0, 2)
    assert list(P.columns) == ["AAPL", "MSFT"]
    assert Q.shape == (0, 1)
    assert list(Q.columns) == ["views"]


def test_sentiment_to_black_litterman_views_validates_input_data():
    namespace = _example_namespace()
    helper = namespace["sentiment_to_black_litterman_views"]

    with pytest.raises(ValueError, match="sentiment_score"):
        helper(pd.DataFrame(index=["AAPL"]), ["AAPL"])

    sentiment = pd.DataFrame({"sentiment_score": [0.1, 0.2]}, index=["AAPL", "AAPL"])
    with pytest.raises(ValueError, match="duplicate tickers"):
        helper(sentiment, ["AAPL"])
