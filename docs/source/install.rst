#######
Install
#######

.. raw:: html

    <a href='https://ko-fi.com/B0B833SXD' target='_blank'><img height='36'style='border:0px;height:36px;' src='https://cdn.ko-fi.com/cdn/kofi1.png?v=2' border='0' alt='Buy Me a Coffee at ko-fi.com' /></a>

Mac OS X, Windows, and Linux
============================

Riskfolio-lib only supports Python 3.7+ on OS X, Windows, and Linux. I recommend
using pip for installation.

1. It is highly recommendable that you must have installed a scientific Python distribution like `anaconda <https://www.anaconda.com/products/individual>` or `winpython <https://winpython.github.io>` (Windows only).

2. If you don't have installed cvxpy, you must follow cvxpy installation `instructions <https://www.cvxpy.org/install/index.html>`_ before installing Riskfolio-Lib.

3. If you still have problems installing cvxpy, you can download cvxpy wheel from this `link <link https://www.lfd.uci.edu/~gohlke/pythonlibs/#cvxpy>` and install using pip.

  ::

      pip install path/cvxpy‑version.whl


4. Install ``Riskfolio-lib``.

  ::

      pip install riskfolio-lib
  

Dependencies
============

Riskfolio-Lib has the following dependencies:

* numpy>=1.17.0
* scipy>=1.1.0
* pandas>=1.0.0
* matplotlib>=3.3.0
* cvxpy>=1.0.25
* scikit-learn>=0.22.0
* statsmodels>=0.10.1
* arch>=4.15
* xlsxwriter>=1.3.7
