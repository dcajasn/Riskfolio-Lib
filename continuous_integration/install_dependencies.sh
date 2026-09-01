#!/bin/bash
# This script is meant to be called by the "Install" step defined in
# build.yml. The behavior of the script is controlled by environment
# variables defined in the build.yml in .github/workflows/.

set -e
uv venv
if [[ "$RUNNER_OS" == "Windows" ]]; then
  . .venv/Scripts/activate
else
  . .venv/bin/activate
fi

#if  [[ "$RUNNER_OS" == "Linux" ]]; then
#  sudo apt install openblas
#fi

uv pip install pytest pytest-cov hypothesis "setuptools>65.5.1"

uv pip install scs clarabel osqp

if [[ "$RUNNER_OS" != "macOS" ]]; then
  uv pip install mkl
fi

# The runtime dependencies and the package itself. Without these the test step
# cannot even import riskfolio. The install is editable so that the compiled
# extension is built into ./riskfolio/external, which is the copy that
# `pytest tests/` imports when run from the repository root.
uv pip install -r requirements.txt
uv pip install -e .
