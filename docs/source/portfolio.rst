######################
Portfolio Optimization
######################

.. raw:: html

    <a href='https://ko-fi.com/B0B833SXD' target='_blank'><img height='36'style='border:0px;height:36px;' src='https://cdn.ko-fi.com/cdn/kofi1.png?v=2' border='0' alt='Buy Me a Coffee at ko-fi.com' /></a>

Some Theory
===========

Mean Risk Portfolio Optimization
--------------------------------

Riskfolio-Lib allows to calculate optimum portfolios that results from optimize
one of the following 4 objective functions:

- **Maximum Return Portfolio:**

.. math::
    \begin{aligned}
    &\underset{w}{\max} & & R (w)\\
    &\text{s.t.} & & Aw \geq B\\
    & & &\phi_{i}(w) \leq c_{i} \; \forall \; i \; \in \; [1,13] \\
    & & & R (w) \geq \overline{\mu}
    \end{aligned}

- **Minimum Risk Portfolio:**

.. math::
    \begin{aligned}
    &\underset{w}{\min} & & \phi_{k}(w)\\
    &\text{s.t.} & & Aw \geq B\\
    & & &\phi_{i}(w) \leq c_{i} \; \forall \; i \; \in \; [1,13] \\
    & & & R (w) \geq \overline{\mu}
    \end{aligned}


- **Maximum Risk Adjusted Return Ratio Portfolio:**

.. math::
    \begin{aligned}
    &\underset{w}{\max} & & \frac{R (w) - r_{f}}{\phi_{k}(w)}\\
    &\text{s.t.} & & Aw \geq B\\
    & & &\phi_{i}(w) \leq c_{i} \; \forall \; i \; \in \; [1,13] \\
    & & & R (w) \geq \overline{\mu}
    \end{aligned}


- **Maximum Utility Portfolio:**

.. math::
    \begin{aligned}
    &\underset{w}{\max} & & R (w) - \lambda \phi_{k}(w)\\
    &\text{s.t.} & & Aw \geq B\\
    & & &\phi_{i}(w) \leq c_{i} \; \forall \; i \; \in \; [1,13] \\
    & & & R (w) \geq \overline{\mu}
    \end{aligned}


Where:

:math:`R (w)` is the return function, posible values are:

    - :math:`\mu w`: arithmetic return.
    - :math:`\mu w - 0.5 w^{\tau} \Sigma w`: approximate logarithmic return :cite:`a-Thorp`.
    - :math:`\frac{1}{T} \sum^{T}_{i=1} \ln (1+ r_{i} w)`: exact logarithmic return :cite:`a-Cajas2`.

:math:`w`: is the vector of weights of the optimum portfolio.

:math:`\mu`: is the vector of expected returns.

:math:`\Sigma`: is the covariance matrix of assets returns.

:math:`r`: is the matrix of assets returns.

:math:`Aw \geq B`: is a set of linear constraints.

:math:`\phi_{i}(w)`: are 18 available risk measures. The available risk
measures are:

- Standard Deviation :cite:`a-Markowitz`.
- Mean Absolute Deviation :cite:`a-Konno`.
- Gini Mean Difference :cite:`a-Yitzhaki1`, :cite:`a-Cajas3`.
- Conditional Value at Risk Range :cite:`a-Cajas3`.
- Tail Gini Range :cite:`a-Cajas3`.
- Semi Standard Deviation :cite:`a-Mansini3`.
- First Lower Partial Moment (Omega Ratio) :cite:`a-Mansini1`.
- Second Lower Partial Moment (Sortino Ratio) :cite:`a-Mansini1`.
- Conditional Value at Risk :cite:`a-Rockafellar`.
- Tail Gini :cite:`a-Ogryczak2002`, :cite:`a-Cajas3`.
- Entropic Value at Risk :cite:`a-Ahmadi2012`, :cite:`a-Ahmadi2017`, :cite:`a-Cajas1`.
- Worst Realization (Minimax) :cite:`a-Mansini2`.
- Maximum Drawdown of uncompounded cumulative returns (Calmar Ratio) :cite:`a-Uryasev1`.
- Average Drawdown of uncompounded cumulative returns :cite:`a-Uryasev1`.
- Conditional Drawdown at Risk of uncompounded cumulative returns :cite:`a-Uryasev1`.
- Entropic Drawdown at Risk of uncompounded cumulative returns :cite:`a-Cajas1`.
- Ulcer Index of uncompounded cumulative returns :cite:`a-martin1989`.

:math:`c_{i}`: are maximum values on each risk measure.

:math:`r_{f}`: is the risk free rate. When the risk measure is the first or
second lower partial moment, :math:`r_{f}` is the minimum acceptable return
:math:`\text{MAR}`.

:math:`\lambda`: is the risk aversion coefficient of the investor.


Risk Parity Portfolio Optimization
----------------------------------

Riskfolio-Lib allows to calculate optimum portfolios that results from optimize
the general vanilla risk parity model :cite:`a-Roncalli` :cite:`a-RichardRoncalli`:

.. math::
    \begin{aligned}
    &\underset{w}{\min} & & \phi(w) \\
    &\text{s.t.} & & b \log(w) \geq c \\
    & & & \mu w \geq \overline{\mu} \\
    & & & Aw \geq B \\
    & & & w \geq 0 \\
    \end{aligned}


Where:
    
:math:`w`: is the vector of weights of the optimum portfolio.

:math:`\mu`: is the vector of expected returns.

:math:`b`: is a vector of risk contribution constraints.

:math:`Aw \geq B`: is a set of linear constraints.

:math:`\phi(w)`: are 10 available risk measures. The available risk
measures are:

- Standard Deviation :cite:`a-Markowitz`.
- Mean Absolute Deviation :cite:`a-Konno`.
- Gini Mean Difference :cite:`a-Yitzhaki1`, :cite:`a-Cajas3`.
- Conditional Value at Risk Range :cite:`a-Cajas3`.
- Tail Gini Range :cite:`a-Cajas3`.
- Semi Standard Deviation :cite:`a-Mansini3`.
- First Lower Partial Moment (Omega Ratio) :cite:`a-Mansini1`.
- Second Lower Partial Moment (Sortino Ratio) :cite:`a-Mansini1`.
- Conditional Value at Risk :cite:`a-Rockafellar`.
- Tail Gini :cite:`a-Ogryczak2002`, :cite:`a-Cajas3`.
- Entropic Value at Risk :cite:`a-Ahmadi2012`, :cite:`a-Ahmadi2017`, :cite:`a-Cajas1`.
- Conditional Drawdown at Risk of uncompounded cumulative returns :cite:`a-Uryasev1`.
- Entropic Drawdown at Risk of uncompounded cumulative returns :cite:`a-Cajas1`.
- Ulcer Index of uncompounded cumulative returns :cite:`a-martin1989`.

:math:`c`: is an arbitrary constant.


Relaxed Risk Parity Portfolio Optimization
------------------------------------------

Riskfolio-Lib allows to calculate optimum portfolios that results from optimize
the relaxed risk parity model :cite:`a-GambetaKwon`:

.. math::
    \begin{aligned}
    &\underset{w}{\min} & & \psi - \gamma & \\
    &\text{s.t.} & & \zeta = \Sigma w \\
    & & & w^{T} \Sigma w \leq N \left ( \psi^{2} - \rho^{2} \right ) & \\
    & & & w_{i} \zeta_{i} \geq \gamma^{2} & \forall i=1 , \ldots , N \\
    & & & \lambda x^{T} \Theta x \leq \rho^{2} & \\
    & & & \mu w \geq \overline{\mu} & \\
    & & & Aw \geq B & \\
    & & & \sum^{N}_{i=1} w_{i} = 1 & \\
    & & & \psi, \gamma, \rho, w  \geq 0 & \\
    \end{aligned}


Where:
    
:math:`w`: is the vector of weights of the optimum portfolio.

:math:`\mu`: is the vector of expected returns.

:math:`\Sigma`: is the covariance matrix of assets returns.

:math:`\psi`: is the average risk of the portfolio.

:math:`\gamma`: is the lower bound of each asset risk constribution.

:math:`\zeta_{i}`: is the marginal risk of asset :math:`i`.

:math:`\rho`: is a regularization variable.

:math:`\lambda`: is a penalty parameter of :math:`\rho`.

:math:`\Theta = \text{diag}(\Sigma)`

:math:`Aw \geq B`: is a set of linear constraints.


Worst Case Mean Variance Portfolio Optimization
-----------------------------------------------

Riskfolio-Lib allows to calculate worst case mean variance optimum portfolios
:cite:`a-Lobo` :cite:`a-fabozzi2007robust` :cite:`a-Tutuncu` :cite:`a-Palomar` 
that results from optimize one of the following 4 objective functions:
       

- **Worst Case Maximum Return Portfolio:**

.. math::
    \begin{aligned}
    &\underset{w}{\max} & & \underset{\mu \, \in \, U_{\mu}}{\min} \mu w\\
    &\text{s.t.} & & Aw \geq B\\
    \end{aligned}


- **Worst Case Minimum Risk Portfolio:**

.. math::
    \begin{aligned}
    &\underset{w}{\max} & & \underset{\Sigma \, \in \, U_{\Sigma}}{\max} w^{T} \Sigma w\\
    &\text{s.t.} & & Aw \geq B\\
    \end{aligned}


- **Worst Case Maximum Risk Adjusted Return Ratio Portfolio:**

.. math::
    \begin{aligned}
    &\underset{w}{\max} & & \cfrac{\underset{\mu \, \in \, U_{\mu}}{\min} \mu w - r_{f}}
    {\underset{\Sigma \, \in \, U_{\Sigma}}{\max} \sqrt{w^{T} \Sigma w}}\\
    &\text{s.t.} & & Aw \geq B\\
    \end{aligned}


- **Worst Case Maximum Utility Portfolio:**

.. math::
    \begin{aligned}
    &\underset{w}{\max} & & \underset{\mu \, \in \, U_{\mu}}{\min} \mu w
    - \underset{\Sigma \, \in \, U_{\Sigma}}{\max}  \lambda w^{T} \Sigma w\\
    &\text{s.t.} & & Aw \geq B\\
    \end{aligned}


Where:

:math:`w` are the weights of the portfolio.

:math:`\mu`: is the vector of expected returns.

:math:`\Sigma` is the covariance matrix.

:math:`U_{\mu}` is the uncertainty set of the mean vector. The uncertainty sets can be:

.. math::
    \begin{aligned}
    U^{box}_{\mu} & = \left \{ \mu \, | \, | \mu - \hat{\mu} | \leq \delta \right \} \\
    U^{ellip}_{\mu} & = \left \{ \mu \, | \left ( \mu - \hat{\mu} \right ) \Sigma^{-1}_{\mu} \left ( \mu - \hat{\mu} \right )^{T} \leq k^{2}_{\mu} \right \} \\
    \end{aligned}


:math:`U_{\Sigma}` is the uncertainty set of the covariance matrix. The uncertainty sets can be:

.. math::
    \begin{aligned}
    U^{box}_{\Sigma} & = \left \{ \Sigma \, | \, \Sigma_{lower} \leq \Sigma \leq \Sigma_{upper} \, , \, \Sigma \succeq 0 \right \} \\
    U^{ellip}_{\Sigma} & = \left \{ \Sigma \, | \left ( \text{vec}(\Sigma) - \text{vec}(\hat{\Sigma}) \right ) \Sigma^{-1}_{\Sigma} \left ( \text{vec}(\Sigma) - \text{vec}(\hat{\Sigma}) \right )^{T} \leq k^{2}_{\Sigma} \, , \, \Sigma \succeq 0 \right \}  \\
    \end{aligned}


:math:`Aw \geq B`: is a set of linear constraints.

:math:`r_{f}`: is the risk free rate. 

:math:`\lambda`: is the risk aversion coefficient of the investor.


Ordered Weighted Averaging (OWA) Portfolio
------------------------------------------

- **Minimum Risk Portfolio:**

.. math::
    \begin{aligned}
    &\underset{w}{\min} & & \sum^{T}_{i=0} v_{[i]}y_{[i]} \\
    &\text{s.t.} & & Aw \geq B\\
    & & & y = rw \\
    & & & R (w) \geq \overline{\mu} \\
    \end{aligned}

- **Maximum Risk Adjusted Return Ratio Portfolio:**

.. math::
    \begin{aligned}
    &\underset{w}{\max} & & \frac{R (w) - r_{f}}{\sum^{T}_{i=0} v_{[i]}y_{[i]}}\\
    &\text{s.t.} & & Aw \geq B\\
    & & & y = rw \\
    & & & R (w) \geq \overline{\mu} \\
    \end{aligned}

- **Maximum Utility Portfolio:**

.. math::
    \begin{aligned}
    &\underset{w}{\max} & & R (w) - \lambda \left ( \sum^{T}_{i=0} v_{[i]}y_{[i]} \right) \\
    &\text{s.t.} & & Aw \geq B\\
    & & & y = rw \\
    & & & R (w) \geq \overline{\mu} \\
    \end{aligned}

Where:

:math:`w` are the weights of the portfolio.

:math:`v` are the weights of the owa operator.

:math:`\mu`: is the vector of expected returns.

:math:`X_{[i]}`: is the element of order :math:`i` of vector :math:`X`.


Module Methods
==============

.. automodule:: Portfolio
   :members:
   :private-members:


Bibliography
============

.. bibliography:: biblio.bib
   :style: unsrt
   :labelprefix: A
   :keyprefix: a-