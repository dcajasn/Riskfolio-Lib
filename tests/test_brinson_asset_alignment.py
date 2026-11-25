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


def test_weight_order_does_not_change_attribution(inputs):
    prices, w, wb, classes = inputs
    expected, _ = attribution(prices, w, wb, classes)
    actual, _ = attribution(prices, w.iloc[::-1], wb, classes)
    pd.testing.assert_frame_equal(actual, expected)
    assert actual.loc["Total Excess Return", "Total"] == pytest.approx(-0.001)


def test_class_order_does_not_change_attribution(inputs):
    prices, w, wb, classes = inputs
    expected, _ = attribution(prices, w, wb, classes)
    actual, _ = attribution(prices, w, wb, classes.iloc[::-1])
    pd.testing.assert_frame_equal(actual, expected)


def test_benchmark_row_vector_is_aligned(inputs):
    prices, w, wb, classes = inputs
    wb = pd.Series([0.3, 0.2, 0.35, 0.15], index=prices.columns)
    expected, _ = attribution(prices, w, wb, classes)
    actual, _ = attribution(prices, w, wb.iloc[::-1].to_frame().T, classes)
    pd.testing.assert_frame_equal(actual, expected)


def test_mismatched_weight_labels_rejected(inputs):
    prices, w, wb, classes = inputs
    w.index = ["A", "B", "C", "UNKNOWN"]
    with pytest.raises(ValueError, match="asset labels"):
        attribution(prices, w, wb, classes)


def test_unlabeled_weights_remain_positional(inputs):
    prices, w, wb, classes = inputs
    expected, _ = attribution(prices, w, wb, classes)
    actual, _ = attribution(prices, w.reset_index(drop=True), wb, classes)
    pd.testing.assert_frame_equal(actual, expected)
