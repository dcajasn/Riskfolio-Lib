#####################
Parameters Estimation
#####################

.. raw:: html

    <a href='https://ko-fi.com/B0B833SXD' target='_blank'><img height='36'style='border:0px;height:36px;' src='https://cdn.ko-fi.com/cdn/kofi1.png?v=2' border='0' alt='Buy Me a Coffee at ko-fi.com' /></a>

This module containts functions that allows us to estimate the vector of means 
and covariance matrix using several methods:

- Historical estimates.
- Estimates using exponencial weighted moving averages (EWMA).
- Robust estimates of the covariance matrix like shrinkage or oracle method.
- The Black Litterman model.
- Factors models to estimate the vector of means and covariance matrix.
- Bootstrapping methods to estimate the input parameters of the uncertainty sets on mean vector and covariance matrix for worst case optimization models.


Module Functions
================

.. automodule:: ParamsEstimation
   :members:
   :private-members:
   :member-order: bysource


Bibliography
============
.. bibliography:: biblio.bib
   :style: unsrt
   :labelprefix: B
   :keyprefix: b-