import numpy as np
import pandas as pd
import pytest

import riskfolio.src.ParamsEstimation as pe


@pytest.mark.parametrize("target", ["b1", "b2", "b3"])
def test_bop_mean_vector_target_normalization(target):
    # Twenty observations: mean (1, 2, 3, 4, 5)/100, covariance I/475.
    # With the existing sample-covariance convention, eqs. (6)-(8) of
    # Bodnar, Okhrin and Parolya give alpha = 17/57 and the vector below.
    # All three targets are nonzero multiples of the same direction;
    # rescaling that direction must not change beta * target.
    mean = np.arange(1, 6) / 100
    deviations = np.vstack([np.eye(5), -np.eye(5)]) / 10
    returns = pd.DataFrame(
        mean + np.tile(deviations, (2, 1)), columns=list("ABCDE")
    )
    expected = np.array([[137, 154, 171, 188, 205]]) / 5700

    result = pe.mean_vector(returns, method="BOP", target=target)

    assert result.columns.tolist() == returns.columns.tolist()
    np.testing.assert_allclose(result.to_numpy(), expected, rtol=1e-12, atol=1e-14)
