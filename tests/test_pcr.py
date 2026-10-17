import numpy as np
import pandas as pd

import riskfolio as rp


def test_pcr_with_all_components_matches_ols():
    rng = np.random.default_rng(0)
    X = pd.DataFrame(rng.normal(0, 0.01, (30, 3)), columns=["F1", "F2", "F3"])
    y = pd.Series(
        0.001 + X.to_numpy() @ np.array([0.5, -1.0, 2.0]) + rng.normal(0, 0.001, 30)
    )

    beta = rp.PCR(X, y, n_components=3)

    ols, *_ = np.linalg.lstsq(
        np.column_stack([np.ones(len(X)), X.to_numpy()]), y.to_numpy(), rcond=None
    )
    np.testing.assert_allclose(beta.ravel(), ols, rtol=1e-10, atol=1e-12)
