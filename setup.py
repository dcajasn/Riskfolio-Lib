# Copyright (C) 2020-2022 Dany Cajas

import os
import numpy as np

from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup

DESCRIPTION = "Portfolio Optimization and Quantitative Strategic Asset Allocation in Python"

with open("README.md", encoding='UTF-8') as fh:
    LONG_DESCRIPTION = fh.read()

DISTNAME = 'Riskfolio-Lib'
MAINTAINER = 'Dany Cajas'
MAINTAINER_EMAIL = 'dany.cajas.n@uni.pe'
URL = 'https://github.com/dcajasn/Riskfolio-Lib'
LICENSE = 'BSD (3-clause)'
KEYWORDS = 'finance, portfolio, optimization, quant, asset, allocation, investing'
DOWNLOAD_URL = 'https://github.com/dcajasn/Riskfolio-Lib.git'
VERSION = '4.0.0'
PYTHON_REQUIRES = ">=3.7"

INSTALL_REQUIRES = [
    'numpy>=1.17.0',
    'scipy>=1.0.1',
    'pandas>=1.0.0',
    'matplotlib>=3.3.0',
    'cvxpy>=1.0.25',
    'scikit-learn>=0.22.0',
    'statsmodels>=0.10.1',
    'arch>=4.15',
    'xlsxwriter>=1.3.7',
    'networkx>=2.5.1',
    'astropy>=4.3.1',
    'pybind11>=2.10.1',
]

PACKAGES = [
    'riskfolio',
    'riskfolio.src',
    'riskfolio.external',
]

CLASSIFIERS = [
    'Intended Audience :: Financial and Insurance Industry',
    'Intended Audience :: Science/Research',
    'Programming Language :: Python :: 3.7',
    'Programming Language :: Python :: 3.8',
    'Programming Language :: Python :: 3.9',
    'License :: OSI Approved :: BSD License',
    'Topic :: Office/Business :: Financial :: Investment',
    'Topic :: Office/Business :: Financial',
    'Topic :: Scientific/Engineering :: Mathematics',
    'Operating System :: Microsoft',
    'Operating System :: Unix',
    'Operating System :: MacOS'
]


if __name__ == "__main__":

    from setuptools import Extension, setup, find_packages
    import sys

    if sys.version_info[:2] < (3, int(PYTHON_REQUIRES[-1])):
        raise RuntimeError("Riskfolio-Lib requires python " + PYTHON_REQUIRES)

    # Obtain the numpy include directory.  This logic works across numpy versions.
    try:
        numpy_include = np.get_include()
    except AttributeError:
        numpy_include = np.get_numpy_include()

    WIN = sys.platform.startswith("win32")
    armadillo_path = os.path.abspath(os.path.join('.', 'lib', 'armadillo-11.4.1', 'include'))
    carma_path = os.path.abspath(os.path.join('.', 'lib', 'carma-0.6.6', 'include'))
    lpack_path = os.path.abspath(os.path.join('.', 'lib', 'armadillo-11.4.1', 'lib_win64'))
    external_path = os.path.abspath(os.path.join('.', 'riskfolio', 'external'))

    if WIN:
        external_module = Pybind11Extension('riskfolio.external.cppfunctions',
            sources=[os.path.join('riskfolio', 'external', 'cpp_functions_bindings.cpp')],
            include_dirs = [numpy_include, armadillo_path, carma_path, external_path],
            define_macros = [('VERSION_INFO', VERSION)],
            )
    else:
        external_module = Pybind11Extension('riskfolio.external.cppfunctions',
            sources=[os.path.join('riskfolio', 'external', 'cpp_functions_bindings.cpp')],
            include_dirs = [numpy_include, armadillo_path, carma_path, external_path],
            extra_compile_args = ['-DARMA_DONT_USE_WRAPPER'],
            define_macros = [('VERSION_INFO', VERSION)],
            )

    setup(
        name=DISTNAME,
        author=MAINTAINER,
        author_email=MAINTAINER_EMAIL,
        maintainer=MAINTAINER,
        maintainer_email=MAINTAINER_EMAIL,
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        long_description_content_type='text/markdown; charset=UTF-8; variant=GFM',
        license=LICENSE,
        keywords=KEYWORDS,
        url=URL,
        version=VERSION,
        download_url=DOWNLOAD_URL,
        python_requires=PYTHON_REQUIRES,
        install_requires=INSTALL_REQUIRES,
        packages=PACKAGES,
        classifiers=CLASSIFIERS,
        project_urls={"Documentation": "https://riskfolio-lib.readthedocs.io/en/latest/",
                      "Issues": "https://github.com/dcajasn/Riskfolio-Lib/issues",
                      "Personal website": "http://financioneroncios.wordpress.com",
                      },
        cmdclass={"build_ext": build_ext},
        ext_modules=[external_module],
    )