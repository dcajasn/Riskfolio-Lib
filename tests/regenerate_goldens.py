"""Rewrite the recorded weight files used by tests/test_portfolio.py.

The recorded files are a regression baseline, not a specification, so only run
this when a change to the optimization models is deliberate, and read the diff
before committing it.

    python tests/regenerate_goldens.py             # rewrite all of them
    python tests/regenerate_goldens.py HC_HRP.csv  # rewrite one

It imports the same functions the tests call, so a recorded file cannot drift
away from the code that produced it.
"""

# Copyright (c) 2020-2026, Dany Cajas
# All rights reserved.
# This work is licensed under BSD 3-Clause "New" or "Revised" License.
# License available at https://github.com/dcajasn/Riskfolio-Lib/blob/master/LICENSE.txt

import os
import sys

sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))

from test_portfolio import GOLDENS, check_invariants, resource  # noqa: E402


def main(argv):
    names = argv[1:] or sorted(GOLDENS)
    unknown = [n for n in names if n not in GOLDENS]
    if unknown:
        raise SystemExit(f"unknown file(s): {unknown}\nknown: {sorted(GOLDENS)}")

    for name in names:
        w = GOLDENS[name]()
        check_invariants(name, w)
        path = resource(name)
        w.to_csv(path)
        print(f"wrote {path}")


if __name__ == "__main__":
    main(sys.argv)
