""""""  #
"""
Copyright (c) 2020-2026, Dany Cajas
All rights reserved.
This work is licensed under BSD 3-Clause "New" or "Revised" License.
License available at https://github.com/dcajasn/Riskfolio-Lib/blob/master/LICENSE.txt
"""

import doctest
import importlib
import pkgutil

import matplotlib
import pytest

import riskfolio

matplotlib.use("Agg")


def _modules():
    names = [riskfolio.__name__]
    for info in pkgutil.walk_packages(riskfolio.__path__, riskfolio.__name__ + "."):
        names.append(info.name)
    return names


def _doctests():
    finder = doctest.DocTestFinder()
    tests = []
    for name in _modules():
        try:
            module = importlib.import_module(name)
        except ImportError:  # optional compiled extensions
            continue
        for test in finder.find(module, name):
            if test.examples:
                tests.append(test)
    return sorted(tests, key=lambda t: t.name)


DOCTESTS = _doctests()


def test_doctests_were_collected():
    # Guards against the suite silently passing because nothing was found.
    assert len(DOCTESTS) > 0


@pytest.mark.parametrize("test", DOCTESTS, ids=lambda t: t.name)
def test_doctest(test):
    runner = doctest.DocTestRunner(optionflags=doctest.NORMALIZE_WHITESPACE)
    output = []
    runner.run(test, out=output.append)
    assert runner.failures == 0, "".join(output)
