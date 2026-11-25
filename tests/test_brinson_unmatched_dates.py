import numpy as np
import pandas as pd
import pytest

from riskfolio.src.RiskFunctions import BrinsonAttribution


@pytest.fixture
def inputs():
    prices = pd.DataFrame(
        [[100., 100., 100., 100.], [110., 95., 103., 108.]],
        index=pd.to_datetime(["2025-01-03", "2025-01-06"]),
        columns=["A", "B", "C", "D"],
    )
    w = pd.Series([0.4, 0.3, 0.2, 0.1], index=prices.columns)
    wb = pd.Series([0.25] * 4, index=prices.columns)
    classes = pd.DataFrame(
        {"Assets": prices.columns, "Sector": ["Tech", "Tech", "Health", "Health"]}
    )
    return prices, w, wb, classes


def attribution(prices, w, wb, classes, **kwargs):
    return BrinsonAttribution(
        prices, w, wb, kwargs.pop("start", "2025-01-03"),
        kwargs.pop("end", "2025-01-06"), classes, "Sector", **kwargs
    )


@pytest.mark.parametrize("missing", ["start", "end"])
def test_exact_date_must_exist(inputs, missing):
    prices, w, wb, classes = inputs
    with pytest.raises(ValueError, match="date cannot be matched"):
        attribution(prices, w, wb, classes, method=None, **{missing: "2025-01-04"})


def test_nearest_and_exact_valid_dates(inputs):
    prices, w, wb, classes = inputs
    exact, dates = attribution(prices, w, wb, classes, method=None)
    nearest, actual_dates = attribution(prices, w, wb, classes, start="2025-01-04")
    assert dates == actual_dates == ("2025-01-03", "2025-01-06")
    pd.testing.assert_frame_equal(nearest, exact)
