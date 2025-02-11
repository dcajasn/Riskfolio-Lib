#####################
Parameters Estimation
#####################

.. raw:: html

    <a href="https://www.kqzyfj.com/click-101359873-15150084?url=https%3A%2F%2Flink.springer.com%2Fbook%2F9783031843037" target="_blank">
        <button style="padding:10px 20px; font-size:16px; background-color: #FFA500; color:white; border:none; border-radius:5px; cursor:pointer;">
            Buy Advanced Portfolio Optimization Book on Springer
        </button>
    </a>
    <br>
    <br>

.. image:: https://img.shields.io/static/v1?label=Sponsor&message=%E2%9D%A4&logo=GitHub&color=%23fe8e86
 :target: https://github.com/sponsors/dcajasn

.. raw:: html
   
    <br>
   
.. raw:: html

    <a href='https://ko-fi.com/B0B833SXD' target='_blank'><img height='36'style='border:0px;height:36px;' src='https://cdn.ko-fi.com/cdn/kofi1.png?v=2' border='0' alt='Buy Me a Coffee at ko-fi.com' /></a>


This module has functions that allows us to estimate the vector of means, covariance matrix and cokurtosis square matrix using several methods:

- Historical estimates.
- Estimates using exponencial weighted moving averages (EWMA).
- Robust estimates of the covariance matrix like Ledoit and Wolf, Oracle, Shrinkage and Graphical Lasso, j-LoGo :cite:`b-jLogo`, Gerber statistic :cite:`b-Gerber2021` and Denoise :cite:`b-MLforAM` estimators.
- Factors models to estimate the vector of means and covariance matrix.
- The Black Litterman model that allows to incorporate analyst's views on returns in estimates of vector of means and covariance matrix :cite:`b-BlackLitterman` :cite:`b-Black1`.
- The Augmented Black Litterman model that allows to incorporate analyst's views on risk factors in estimates of vector of means and covariance matrix :cite:`b-WCheung`.
- The Black Litterman Bayesian model that allows to incorporate analyst's views on risk factors in estimates of vector of means and covariance matrix :cite:`b-BLB`.
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
