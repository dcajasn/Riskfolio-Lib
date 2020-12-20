################
Portfolio Models
################

Some Theory
===========

Mean Risk Portfolios
--------------------

Riskfolio-Lib allows to calculate optimum portfolios that results from optimize
one of the following 4 objective functions:

- **Maximum Return Portfolio:**

.. math::

    \begin{align}
    &\underset{w}{\max} & & \mu w\\
    &\text{s.t.} & & Aw \geq B\\
    & & &\phi_{i}(w) \leq c_{i} \; \forall \; i \; \in \; [1,10] \\
    \end{align}

- **Minimum Risk Portfolio:**

.. math::

    \begin{align}
    &\underset{w}{\min} & & \phi_{k}(w)\\
    &\text{s.t.} & & Aw \geq B\\
    & & &\phi_{i}(w) \leq c_{i} \; \forall \; i \; \in \; [1,10] \\
    \end{align}

- **Maximum Risk Adjusted Return Ratio Portfolio:**

.. math::

    \begin{align}
    &\underset{w}{\max} & & \frac{\mu w - r_{f}}{\phi_{k}(w)}\\
    &\text{s.t.} & & Aw \geq B\\
    & & &\phi_{i}(w) \leq c_{i} \; \forall \; i \; \in \; [1,10] \\
    \end{align}

- **Maximum Utility Portfolio:**

.. math::

    \begin{align}
    &\underset{w}{\max} & & \mu w - \lambda \phi_{k}(w)\\
    &\text{s.t.} & & Aw \geq B\\
    & & &\phi_{i}(w) \leq c_{i} \; \forall \; i \; \in \; [1,10] \\
    \end{align}

Where:
    
:math:`\mu`: is the vector of expected returns.

:math:`w`: is the vector of weights of the optimum portfolio.

:math:`Aw \geq B`: is a set of linear constraints.

:math:`\phi_{i}(w)`: are 10 available risk measures. The available risk
measures are:

- Standard Deviation :cite:`a-Markowitz`.
- Mean Absolute Deviation :cite:`a-Konno`.
- Semi Standard Deviation :cite:`a-Mansini3`.
- First Lower Partial Moment (Omega Ratio) :cite:`a-Mansini1`.
- Second Lower Partial Moment (Sortino Ratio) :cite:`a-Mansini1`.
- Conditional Value at Risk :cite:`a-Rockafellar`.
- Worst Realization (Minimax) :cite:`a-Mansini2`
- Maximum Drawdown of uncompounded returns (Calmar Ratio) :cite:`a-Uryasev1`.
- Average Drawdown of uncompounded returns :cite:`a-Uryasev1`.
- Conditional Drawdown at Risk of uncompounded returns :cite:`a-Uryasev1`.

:math:`c_{i}`: are maximum values on each risk measure.

:math:`r_{f}`: is the risk free rate. When the risk measure is the first or
second lower partial moment, :math:`r_{f}` is the minimum acceptable return
:math:`\text{MAR}`.

:math:`\lambda`: is the risk aversion coefficient of the investor.


Risk Parity Portfolios
----------------------

Riskfolio-Lib allows to calculate optimum portfolios that results from optimize
the general risk parity model :cite:`a-Roncalli`:

.. math::

    \begin{align}
    &\underset{w}{\min} & & \phi(w)\\
    &\text{s.t.} & & b \log(w) \geq c\\
    & & & w \geq 0 \\
    \end{align}

Where:
    
:math:`w`: is the vector of weights of the optimum portfolio.

:math:`b`: is a vector of risk contribution constraints.

:math:`\phi(w)`: are 7 available risk measures. The available risk
measures are:

- Standard Deviation :cite:`a-Markowitz`.
- Mean Absolute Deviation :cite:`a-Konno`.
- Semi Standard Deviation :cite:`a-Mansini3`.
- First Lower Partial Moment (Omega Ratio) :cite:`a-Mansini1`.
- Second Lower Partial Moment (Sortino Ratio) :cite:`a-Mansini1`.
- Conditional Value at Risk :cite:`a-Rockafellar`.
- Conditional Drawdown at Risk of uncompounded returns :cite:`a-Uryasev1`.

:math:`c`: is an arbitrary constant.

Worst Case Mean Variance Portfolios
-----------------------------------

Riskfolio-Lib allows to calculate worst case mean variance optimum portfolios
:cite:`a-Lobo` :cite:`a-fabozzi2007robust` :cite:`a-Tutuncu` :cite:`a-Palomar` 
that results from optimize one of the following 4 objective functions:
       

- **Worst Case Maximum Return Portfolio:**

.. math::

    \begin{align}
    &\underset{w}{\max} & & \underset{\mu \, \in \, U_{\mu}}{\min} \mu w\\
    &\text{s.t.} & & Aw \geq B\\
    \end{align}

- **Worst Case Minimum Risk Portfolio:**

.. math::

    \begin{align}
    &\underset{w}{\max} & & \underset{\Sigma \, \in \, U_{\Sigma}}{\max} w^{T} \Sigma w\\
    &\text{s.t.} & & Aw \geq B\\
    \end{align}
            
- **Worst Case Maximum Risk Adjusted Return Ratio Portfolio:**

.. math::

    \begin{align}
    &\underset{w}{\max} & & \cfrac{\underset{\mu \, \in \, U_{\mu}}{\min} \mu w - r_{f}}
    {\underset{\Sigma \, \in \, U_{\Sigma}}{\max} \sqrt{w^{T} \Sigma w}}\\
    &\text{s.t.} & & Aw \geq B\\
    \end{align}

- **Worst Case Maximum Utility Portfolio:**

.. math::

    \begin{align}
    &\underset{w}{\max} & & \underset{\mu \, \in \, U_{\mu}}{\min} \mu w
    - \underset{\Sigma \, \in \, U_{\Sigma}}{\max}  \lambda w^{T} \Sigma w\\
    &\text{s.t.} & & Aw \geq B\\
    \end{align}

Where:

:math:`w` are the weights of the portfolio.

:math:`\mu`: is the vector of expected returns.

:math:`\Sigma` is the covariance matrix.

:math:`U_{\mu}` is the uncertainty set of the mean vector. The uncertainty sets can be:

.. math::

    \begin{align}
    U^{box}_{\mu} & = \{ \mu \, | \, | \mu - \hat{\mu} | \geq \delta \} \\
    U^{ellip}_{\mu} & = \{ \mu \, | \mu = \hat{\mu} + \kappa \Sigma^{1/2}_{\mu} \} \\
    \end{align}

:math:`U_{\Sigma}` is the uncertainty set of the covariance matrix.

.. math::

    \begin{align}
    U^{box}_{\Sigma} = \{ \Sigma \, | \, \Sigma_{lower} \geq \Sigma \geq \Sigma_{upper} \} \\
    \end{align}

:math:`Aw \geq B`: is a set of linear constraints.

:math:`r_{f}`: is the risk free rate. 

:math:`\lambda`: is the risk aversion coefficient of the investor.


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
