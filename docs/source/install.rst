#######
Install
#######

.. image:: https://img.shields.io/static/v1?label=Sponsor&message=%E2%9D%A4&logo=GitHub&color=%23fe8e86
 :target: https://github.com/sponsors/dcajasn

.. raw:: html
   
   <br>

.. raw:: html

    <a href='https://ko-fi.com/B0B833SXD' target='_blank'><img height='36'style='border:0px;height:36px;' src='https://cdn.ko-fi.com/cdn/kofi1.png?v=2' border='0' alt='Buy Me a Coffee at ko-fi.com' /></a>

Mac OS X, Windows, and Linux
============================

Riskfolio-lib only supports Python 3.7+ on OS X, Windows, and Linux. I recommend
using pip for installation.

1. It is highly recommendable that you must have installed a scientific Python distribution like `anaconda <https://www.anaconda.com/products/individual>`_ or `winpython <https://winpython.github.io>`_ (Windows only).

2. Install ``Pybind11``.

  ::

      pip install pybind11


3. If you don't have installed cvxpy, you must follow `cvxpy <https://www.cvxpy.org/install/index.html>`_ installation instructions before installing Riskfolio-Lib.

4. If you still have problems installing cvxpy, you can download cvxpy wheel from the `Unofficial Windows Binaries for Python Extension Packages <https://www.lfd.uci.edu/~gohlke/pythonlibs/#cvxpy>`_ and install using pip.

  ::

      pip install path/cvxpy‑version.whl

5. Install `Visual Studio Build Tools <https://visualstudio.microsoft.com/es/downloads/>`_ (Only for Windows).

.. image:: images/MVSC1.png
    :width: 100%

.. image:: images/MVSC2.png
    :width: 100%

6. Install ``Riskfolio-lib``.

  ::

      pip install riskfolio-lib


7. To run some examples is necessary to install `yfinance <https://pypi.org/project/yfinance/>`_.

  ::

      pip install yfinance
  

8. To run some examples is necessary to install MOSEK, you must follow `MOSEK <https://docs.mosek.com/9.2/install/installation.html>`_ installation instructions. To get a MOSEK license you must go to `Academic Licenses <https://www.mosek.com/products/academic-licenses/>`_.

 ::

      pip install mosek

Dependencies
============

Riskfolio-Lib has the following dependencies:

* numpy>=1.18.0
* scipy>=1.1.0
* pandas>=1.0.0
* matplotlib>=3.5.0
* cvxpy>=1.3.1
* scikit-learn>=1.0.0
* statsmodels>=0.10.1
* arch>=4.15
* xlsxwriter>=1.3.7
* networkx>=2.5.1
* astropy>=4.3.1 (if there are problems check `astropy installation instructions <https://www.astropy.org>`_)
* pybind11>=2.10.1
