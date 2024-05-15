# Optimal derivative positioning on a single stock

- Can we outperform a benchmark by trading options on it?

- This question may be answered by maximizing $U(w_1f_1(X_T)+...+w_Nf_N(X_T)-We^{X_T})$ subject to the constraints $\mathbb{E}^{\mathbb{Q}}[w_1f_1(X_T)+...+w_Nf_N(X_T)] \leq W$
- $W$ denotes the available capital to be allocated and it is assumed that, denoting by $\Phi$ the Fenchel conjugate of a given distortion $\Psi$,
  - $U:L^{\infty}\rightarrow \mathbb{R}$ is defined by $U(Y) := \min_{Z\in\mathcal{M}}\mathbb{E}[ZY]-\alpha(Z)$ for $Y\in L^{\infty}$, where $\mathcal{M} := \{Z\in L^1_+:\mathbb{E}[Z]=1,\mathbb{E}[(Z-a)^+]\leq \Phi(a), a > 0\}$ and $\alpha(Z) := \mathbb{E}[\theta Z^{\alpha}+(1-\theta)Z^{-\beta}]$
  - $f_1,...,f_N$ are the respective payoff functions of $N$ contingent claims (e.g. options) each with a fixed maturity $T$ (here, 7 days) on the same underlying asset
  - The random variable $X_T$ is the log returns of the underlying asset, and its distribution is assumed to follow the bilateral gamma distribution with parameters $(b_p,c_p,b_n,c_n)$ and $(\tilde{b}_p,\tilde{c}_p,\tilde{b}_n,\tilde{c}_n)$ under the statistical probability $\mathbb{P}$ and the risk neutral probability $\mathbb{Q}$ respectively
  - The parameters $\theta,\alpha,\beta$ may be estimated based on performance of the resulting trading strategy

- **Note on optimal amount $W$ to be invested**: Two possible ways to move forward instead of maximizing net return
  - Specify an alternative (e.g. Treasuries), which then needs to be included in the optimization (1)
  - Use non monotonic MUF -> this may be the only way if, say, we are considering a team that specializes in trading a particular asset class, and the capital allocated to the team is decided only once per year

## Formulation for discipline saddle programming

- To use the dsp extension of CVXPY we need to discretize the problem

- Given a random variable $Y$, $M_0<M_1\in \mathbb{R}$ and $k\in \mathbb{N}$, define
$y_j := M_0+\delta j, \ j\in\{0,1,...,2^k-1\}$
$\delta & := \frac{M_1-M_0}{2^k-1}$
$Y^k & := \sum_{j=1}^{2^k-1}y_j\mathbb{1}_{[y_j,y_{j+1})}(Y)$

- Then, letting $p_Y$ denote the density of $Y$ under $\mathbb{P}$, $\mathbb{P}\left(Y^k=y_j\right)\approx \delta p_Y\left(y_j\right)$

- One is then led to consider the problem
$\max_{\mathbf{w} \in \mathbb{R}^N} \min_{Z^k\in\mathcal{M}} \mathbb{E}[Z^k(w_1f_1(X^k_T)+...+w_Nf_N(X^k_T)-We^{X^k_T})] - \alpha(Z^k)$,
$\text{s.t. } \mathbb{E}^{\mathbb{Q}}[w_1f_1(X_T)+...+w_Nf_N(X_T)] \leq W$

- Next, assuming all strikes are traded for maturity $T$, we consider the problem,
$\max_{\mathbf{q}^T\mathbf{y}=W}
    \min_{\mathbf{z}\in\mathcal{M}} \mathbf{z}^TP_k(\mathbf{y}-We^{\mathbf{x}}) - \mathbf{p}^T(\theta \mathbf{z}^{\alpha}+(1-\theta)\mathbf{z}^{-\beta})$
- In (2), we denote by
  - $P_k$ denote the $2^{k}\times 2^{k}$ diagonal matrix with diagonal elements given by $\delta p_{X}(x_j)$, $j=1,...,2^k$,
  - $\mathbf{p}$ and $\mathbf{q}$ are the probability mass of $X^k_T$ under the statistical and the risk neutral measure respectively
  - $\mathbf{y}$, $\mathbf{x}$ and $\mathbf{z}$ are vectors in $\mathbb{R}^{2^k}$ composed of all possible values of $Y^k$, $X^k$ and $Z^k$
- This formulation seems to fit the functions available in the DSP extension of CVXPY, as it is shown next

## Disciplined Saddle Programming

- The github repository for the dsp extension of CVXPY is available at <https://github.com/cvxgrp/dsp>

- Relevant references on which the repo is based on are:
  - <https://arxiv.org/abs/2301.13427>
  - <https://arxiv.org/abs/2102.01002>

Below is a simple example that computes the solution to problem (2) in dimension $k = 2$

import dsp
import cvxpy as cp
import numpy as np

y = cp.Variable(2)
z = cp.Variable(2)
x = np.array([10,3])
p = np.array([0.3,0.7])
q = np.array([0.5,0.5])
a = 1
Phi = 3
P = np.array([[p[0], 0], [0, p[1]]])
theta = 0.5
alpha = 0.2
beta = 0.3
W = 1

f = dsp.inner(z, P @ (y - cp.multiply(W,cp.exp(x))))
rho = p @ (cp.multiply(theta, cp.power(z,alpha))+cp.multiply(1-theta,cp.power(z,-beta)))
obj = dsp.MinimizeMaximize(rho-f)
constraints = [q @ y == W, p @ z == 1, z >= 0, p @ cp.maximum(z-a,0) <= Phi]

obj = dsp.MinimizeMaximize(f)
prob = dsp.SaddlePointProblem(obj, constraints)
prob.solve()  # solves the problem
