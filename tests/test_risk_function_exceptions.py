"""Solver failures may fall back; cancellation and unexpected errors must not."""

import cvxpy as cp
import numpy as np
import pytest

import riskfolio as rp
import riskfolio.src.RiskFunctions as rk


@pytest.fixture
def returns():
    return np.array([-0.04, -0.02, 0.01, 0.03]).reshape(-1, 1)


@pytest.mark.parametrize("stage", ["evar", "rlvar_dual", "rlvar_primal"])
@pytest.mark.parametrize(
    "error_type", [KeyboardInterrupt, SystemExit, MemoryError, RuntimeError, ValueError]
)
def test_risk_solver_propagates_unexpected_errors(
    monkeypatch, returns, stage, error_type
):
    error = error_type("solver interrupted")
    calls = []

    def solve(problem, **kwargs):
        calls.append(problem)
        if stage == "rlvar_primal" and isinstance(problem.objective, cp.Maximize):
            raise cp.error.SolverError("dual solver failed")
        raise error

    monkeypatch.setattr(cp.Problem, "solve", solve)
    function = rp.EVaR_Hist if stage == "evar" else rp.RLVaR_Hist

    with pytest.raises(error_type) as caught:
        function(returns, alpha=0.5)

    assert caught.value is error
    assert len(calls) == (2 if stage == "rlvar_primal" else 1)


@pytest.mark.parametrize("failure", ["solver_error", "no_value"])
def test_evar_retains_scipy_fallback(monkeypatch, returns, failure):
    solve_calls = []
    minimize_calls = []
    minimize = rk.minimize

    def solve(problem, **kwargs):
        solve_calls.append(kwargs["solver"])
        if failure == "solver_error":
            raise cp.error.SolverError("solver unavailable")

    def record_minimize(*args, **kwargs):
        minimize_calls.append(kwargs["method"])
        return minimize(*args, **kwargs)

    monkeypatch.setattr(cp.Problem, "solve", solve)
    monkeypatch.setattr(rk, "minimize", record_minimize)

    value, z = rp.EVaR_Hist(returns, alpha=0.5)

    assert minimize_calls == ["SLSQP"]
    assert len(solve_calls) == (1 if failure == "solver_error" else 5)
    assert np.isfinite(value)
    assert -returns.mean() <= value <= -returns.min()
    assert z > 0


@pytest.mark.parametrize("failure", ["solver_error", "no_value"])
def test_rlvar_retains_primal_fallback(monkeypatch, returns, failure):
    expected = rp.RLVaR_Hist(returns, alpha=0.5)
    solve = cp.Problem.solve
    calls = []

    def fail_dual(problem, **kwargs):
        calls.append(type(problem.objective))
        if isinstance(problem.objective, cp.Maximize):
            if failure == "solver_error":
                raise cp.error.SolverError("dual solver failed")
            return None
        return solve(problem, **kwargs)

    monkeypatch.setattr(cp.Problem, "solve", fail_dual)

    actual = rp.RLVaR_Hist(returns, alpha=0.5)

    assert calls == [cp.Maximize] * (1 if failure == "solver_error" else 3) + [
        cp.Minimize
    ]
    assert actual == pytest.approx(expected, abs=1e-6)


@pytest.mark.parametrize("failure", ["solver_error", "no_value"])
def test_rlvar_retains_result_when_both_formulations_fail(
    monkeypatch, returns, failure
):
    calls = []

    def solve(problem, **kwargs):
        calls.append(type(problem.objective))
        if failure == "solver_error":
            raise cp.error.SolverError("solver unavailable")

    monkeypatch.setattr(cp.Problem, "solve", solve)

    assert rp.RLVaR_Hist(returns, alpha=0.5) == 0
    count = 1 if failure == "solver_error" else 3
    assert calls == [cp.Maximize] * count + [cp.Minimize] * count


@pytest.mark.parametrize("function", [rp.EVaR_Hist, rp.RLVaR_Hist])
def test_successful_risk_solve_stops_at_first_solver(monkeypatch, returns, function):
    solve = cp.Problem.solve
    calls = []

    def record_solve(problem, **kwargs):
        calls.append(kwargs["solver"])
        return solve(problem, **kwargs)

    monkeypatch.setattr(cp.Problem, "solve", record_solve)

    result = function(returns, alpha=0.5, solver="CLARABEL")
    value = result[0] if isinstance(result, tuple) else result

    assert calls == ["CLARABEL"]
    assert np.isfinite(value)
    assert -returns.mean() <= value <= -returns.min()
