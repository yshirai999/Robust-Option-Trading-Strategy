# Optimal Spot Slides

- The following problem is here considered:
  - $U(w_1f_1(X_T)+...+w_Nf_N(X_T)-\mathbb{E}[w_1f_1(X_T)+...+w_Nf_N(X_T)])$
- It is assumed that, denoting by $\Phi$ the Fenchel conjugate of a given distortion $\Psi$,
  - $U:L^{\infty}\rightarrow \mathbb{R}$ is defined by
    - $U(Y) := \min_{Z\in\mathcal{M}}\mathbb{E}[ZY]-\alpha(Z)$ for $Y\in L^{\infty}$,
    - $\mathcal{M} := \{Z\in L^1_+:\mathbb{E}[Z]=1,\mathbb{E}[(Z-a)^+]\leq \Phi(a), a > 0\}$
    - $\alpha(Z) := \mathbb{E}[\theta Z^{\alpha}+(1-\theta)Z^{-\beta}]$
  - $f_1,...,f_N$ are the respective payoff functions of $N$ contingent claims (e.g. options) each with a fixed maturity $T$ (here, 7 days) on the same underlying asset
  - The random variable $X_T$ is the log returns of the underlying asset, and its distribution is assumed to follow the bilateral gamma distribution with parameters $(b_p,c_p,b_n,c_n)$ and $(\tilde{b}_p,\tilde{c}_p,\tilde{b}_n,\tilde{c}_n)$ under the statistical probability $\mathbb{P}$ and the risk neutral probability $\mathbb{Q}$ respectively
  - The parameters $\theta,\alpha,\beta$ may be estimated based on performance of the resulting trading strategy

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
- This formulation seems to fit the functions available in the DSP extension of CVXPY

## Disciplined Saddle Programming

- The github repository for the dsp extension of CVXPY is available at <https://github.com/cvxgrp/dsp>

- Relevant references on which the repo is based on are:
  - <https://arxiv.org/abs/2301.13427>
  - <https://arxiv.org/abs/2102.01002>


## How to use the repo

- The MATLAB folder contains a matlab m file, which solves the problem considered by running a python script which calls dsp

- Supporte python version is v3.10.14
