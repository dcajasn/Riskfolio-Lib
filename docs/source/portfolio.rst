################
Portfolio Models
################

Some Theory
===========

Riskfolio-Lib allows to calculate optimum portfolios that result from optimize
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
