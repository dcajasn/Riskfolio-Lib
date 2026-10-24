import numpy as np
import pytest

from riskfolio.src.OwaWeights import owa_rg
from riskfolio.src.RiskFunctions import RG


@pytest.mark.parametrize(
    "returns",
    [[0.1], [-0.3], [0.0], [0.1, 0.2], [-0.3, -0.1, 0.2], [0.2, 0.2]],
)
def test_owa_range_matches_peak_to_peak(returns):
    values = np.asarray(returns).reshape(-1, 1)
    weights = owa_rg(len(values))
    assert weights.shape == values.shape
    assert weights.sum() == 0
    result = (weights.T @ np.sort(values, axis=0)).item()
    assert result == pytest.approx(np.ptp(values))
    assert RG(values) == pytest.approx(np.ptp(values))
