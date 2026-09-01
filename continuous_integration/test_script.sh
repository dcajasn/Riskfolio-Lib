#!/bin/bash
# This script is meant to be called by the "Test" step defined in build.yml.
# The behavior of the script is controlled by environment variables defined in
# the build.yml in .github/workflows/.

set -e
if [[ "$RUNNER_OS" == "Windows" ]]; then
  . .venv/Scripts/activate
else
  . .venv/bin/activate
fi

# Coverage is only collected on the single cell that uploads it, so the other
# cells are not slowed down by tracing.
if [[ "$SINGLE_ACTION_CONFIG" == "True" ]]; then
  python -m pytest tests/ -v --cov=riskfolio --cov-report=xml --cov-report=term
else
  python -m pytest tests/ -v
fi
