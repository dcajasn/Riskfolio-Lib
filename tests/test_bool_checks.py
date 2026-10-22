import numpy as np
import pandas as pd

import riskfolio as rp
import riskfolio.src.AuxFunctions as af


def test_is_pos_def_returns_python_bool():
    assert af.is_pos_def(np.eye(3)) is True
    assert af.is_pos_def(np.array([[1.0, 2.0], [2.0, 1.0]])) is False


def test_hrp_constraints_applies_enabled_constraints():
    asset_classes = pd.DataFrame({"Assets": ["A", "B"]})
    constraints = pd.DataFrame(
        {
            "Disabled": [False],
            "Type": ["Assets"],
            "Set": [""],
            "Position": ["A"],
            "Sign": ["<="],
            "Weight": [0.3],
        }
    )
    w_max, _ = rp.hrp_constraints(constraints, asset_classes)
    assert w_max["A"] == 0.3
