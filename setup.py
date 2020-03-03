# Copyright (C) 2020 Dany Cajas

DESCRIPTION = "Riskfolio-Lib: Quantitative Strategic Asset Allocation, easy for you"

with open("README.md", "r") as fh:
    LONG_DESCRIPTION = fh.read()

DISTNAME = 'Riskfolio-Lib'
MAINTAINER = 'Dany Cajas'
MAINTAINER_EMAIL = 'dany.cajas.n@uni.pe'
URL = 'https://riskfolio-lib.readthedocs.io/en/latest/'
LICENSE = 'BSD (3-clause)'
DOWNLOAD_URL = 'https://github.com/dcajasn/Riskfolio-Lib.git'
VERSION = '0.0.2'
PYTHON_REQUIRES = ">=3.6"

INSTALL_REQUIRES = [
    'numpy>=1.17.0',
    'scipy>=1.0.1',
    'pandas>=1.0.0',
    'matplotlib>=3.0.0',
    'cvxpy>=1.0.25',
    'scikit-learn>=0.22.0',
    'statsmodels>=0.10.1',
]


PACKAGES = [
    'riskfolio',
]

CLASSIFIERS = [
    'Intended Audience :: Financial and Insurance Industry',
    'Intended Audience :: Science/Research',
    'Programming Language :: Python :: 3.6',
    'Programming Language :: Python :: 3.7',
    'Programming Language :: Python :: 3.8',
    'License :: OSI Approved :: BSD License',
    'Topic :: Office/Business :: Financial :: Investment',
    'Topic :: Office/Business :: Financial',
    'Operating System :: Microsoft',
    'Operating System :: Unix',
    'Operating System :: MacOS'
]


if __name__ == "__main__":

    from setuptools import setup

    import sys
    if sys.version_info[:2] < (3, 6):
        raise RuntimeError("Riskfolio-Lib requires python >= 3.6.")

    setup(
        name=DISTNAME,
        author=MAINTAINER,
        author_email=MAINTAINER_EMAIL,
        maintainer=MAINTAINER,
        maintainer_email=MAINTAINER_EMAIL,
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        long_description_content_type="text/markdown",
        license=LICENSE,
        url=URL,
        version=VERSION,
        download_url=DOWNLOAD_URL,
        python_requires=PYTHON_REQUIRES,
        install_requires=INSTALL_REQUIRES,
        packages=PACKAGES,
        classifiers=CLASSIFIERS
    )