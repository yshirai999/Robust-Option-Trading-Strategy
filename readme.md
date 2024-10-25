# Portfolio Theory with Monetary Preferences

- The following problem is here considered:
  - $\max_f U(f(X)-\mathbb{E}[f(X)])$
- It is assumed that, denoting by $\Phi$ the Fenchel conjugate of a given distortion $\Psi$,
  - $U:L^{\infty}\rightarrow \mathbb{R}$ is defined by
    - $U(Y) := \min_{Z\in\mathcal{M}}\mathbb{E}[ZY]-\alpha(Z)$ for $Y\in L^{\infty}$,
    - $\mathcal{M} := \{Z\in L^1_+:\mathbb{E}[Z]=1,\mathbb{E}[(Z-a)^+]\leq \Phi(a), a > 0\}$
    - $\alpha(Z) := \mathbb{E}[\theta Z^{\alpha}+(1-\theta)Z^{-\beta}]$
  - The random variable $X$ is the log returns of an underlying risk which can be traded and is liquid
  - The distribution of $X$ is assumed to follow the bilateral gamma distribution with parameters $(b_p,c_p,b_n,c_n)$ under the risk neutral probability $\mathbb{Q}$ 
  - The parameters $\theta,\alpha,\beta$ may be estimated based on performance of the resulting trading strategy

## Idea
- The functional $U$ is the infimum over a set of test measures that are equivalent to a base measure $\mathbb{Q}'$
- Departures from $\mathbb{Q}'$ are penalized
- The time at which $f(X)$ is sold is very short, say a day or so, in such a way that options' moneyness remains unchanged
- Thus, if the underlying options traded have expiration, say, 5 days, then $\mathbb{Q}$ is the risk neutral measure calibrated to options expiring in 5 days
- If the plan is to sell $f(X)$ in 1 day, then $\mathbb{Q}'$ is the risk neutral measure calibrated to options expiring in 4 days
- The distortion of $\mathbb{Q}'$ is justfied because of fluctuations in tomorrow's risk neutral measure compared to today's

## Formulation for discipline saddle programming

- To use the dsp extension of CVXPY we need to discretize the problem

- Given a random variable $Y$, real numbers $M_0\leq M_1$ and an positive integer $k$, define, for $j = 0,1,...,2^k-1$,
  - $y_j := M_0+\delta j$
  - $A_j := [y_j,y_{j+1}]$
  - $\delta := \frac{M_1-M_0}{2^k-1}$
  - $Y^k := \sum_{j}y_j\mathbb{1}_{A_j}(Y)$

- Then, letting $p_Y$ denote the density of $Y$ under $\mathbb{P}$, $\mathbb{P}\left(Y^k=y_j\right)\approx \delta p_Y\left(y_j\right)$

- One is then led to consider the problem
  - $\max_{\mathbf{w} \in \mathbb{R}^N} \min_{Z^k\in\mathcal{M}} \mathbb{E}[Z^k(w_1f_1(X^k_T)+...+w_Nf_N(X^k_T)-We^{X^k_T}) - \mathbb{E}^{\mathbb{Q}}[w_1f_1(X_T)+...+w_Nf_N(X_T)]] + \alpha(Z^k)$

- Next, assuming all strikes are traded for maturity $T$, we consider the problem,
  - $\max_{\mathbf{q}^T\mathbf{y}=W} \min_{\mathbf{z}\in\mathcal{M}} \mathbf{z}^TP_k(\mathbf{y}-We^{\mathbf{x}}) - \mathbf{p}^T(\theta \mathbf{z}^{\alpha}+(1-\theta)\mathbf{z}^{-\beta})$
- where:
  - $P_k$ denotes the $2^{k}\times 2^{k}$ diagonal matrix with diagonal elements given by $\delta p_{X}(x_j)$, $j=1,...,2^k$,
  - $\mathbf{p}$ and $\mathbf{q}$ are the probability mass of $X^k_T$ under the statistical and the risk neutral measure respectively
  - $\mathbf{y}$, $\mathbf{x}$ and $\mathbf{z}$ are vectors in $\mathbb{R}^{2^k}$ composed of all possible values of $Y^k$, $X^k$ and $Z^k$

- The github repository for the dsp extension of CVXPY is available at <https://github.com/cvxgrp/dsp>

- Relevant references on which the repo is based on are:
  - <https://arxiv.org/abs/2301.13427>
  - <https://arxiv.org/abs/2102.01002>
