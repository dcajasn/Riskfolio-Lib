#####################
Constraints Functions
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


This module has functions that help us to create any kind of linear constraint
related to the assets or assets class weights or related to the value of the sensitivity
of the portfolio to a specific risk factor. These functions transform all
constraint to the form :math:`Aw \geq B`.

This module have a function that help us to create relative and absolute views for the Black
Litterman model :cite:`e-BlackLitterman`. This views can consider relationships among assets and asset classes. This
function transform all views to the form :math:`Pw = Q`.

This module also have functions to create constraints based on graph information :cite:`e-Cajas10` :cite:`e-Cajas11`
like the information obtained from networks and dendrograms.

Module Functions
================

.. automodule:: ConstraintsFunctions
   :members:
   :private-members:

Bibliography
============

.. bibliography:: biblio.bib
   :style: unsrt
   :labelprefix: E
   :keyprefix: e-
