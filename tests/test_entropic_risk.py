import numpy as np
import pytest

from riskfolio.src.RiskFunctions import Entropic_RM, _Entropic_RM


@pytest.mark.parametrize("private", [False, True])
@pytest.mark.parametrize("constant", [-1.0, 1.0])
def test_entropic_risk_large_scaled_constant(private, constant):
    z, alpha = 0.001, 0.05
    returns = np.full(4, constant)
    with np.errstate(over="raise", divide="raise", under="ignore"):
        result = (
            _Entropic_RM(z, returns, alpha)
            if private
            else Entropic_RM(returns, z, alpha)
        )
    expected = -constant - z * np.log(alpha)
    assert result == pytest.approx(expected, rel=1e-12)


@pytest.mark.parametrize("private", [False, True])
def test_entropic_risk_large_scaled_mixed_returns(private):
    z, alpha = 0.001, 0.05
    returns = np.array([-1.0, -0.999])
    with np.errstate(over="raise", divide="raise", under="ignore"):
        result = (
            _Entropic_RM(z, returns, alpha)
            if private
            else Entropic_RM(returns, z, alpha)
        )
    expected = 1 + z * (np.log1p(np.exp(-1)) - np.log(2) - np.log(alpha))
    assert result == pytest.approx(expected, rel=1e-12)


@pytest.mark.parametrize("private", [False, True])
def test_entropic_risk_ordinary_values(private):
    returns = np.array([-0.02, 0.01, 0.03])
    z, alpha = 0.5, 0.05
    expected = z * np.log(np.mean(np.exp(-returns / z)) / alpha)
    result = (
        _Entropic_RM(z, returns, alpha)
        if private
        else Entropic_RM(returns, z, alpha)
    )
    assert result == pytest.approx(expected, rel=1e-12)
