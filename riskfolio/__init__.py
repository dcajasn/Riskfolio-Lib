# -*- coding: utf-8 -*-
from pkg_resources import get_distribution, DistributionNotFound

try:
    # Change here if project is renamed and does not equal the package name
    dist_name = "Riskfolio-Lib"
    __version__ = get_distribution(dist_name).version
except DistributionNotFound:
    __version__ = "0.1"
finally:
    del get_distribution, DistributionNotFound
