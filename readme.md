# Optimal Spot Slide

- A desk holding option positions on a single underlying stock with a current price of $S_0$ constructs a spot slide $\zeta(x)$ which represents the logarithm of the aggregate value of all positions were the market price of the underlying instantaneously jumps to the level $S=S_0e^x$.

- Such calculations are now routinely mandated by Central Banks monitoring the health of the banking system. Existing exposures may be altered to match optimal ones, or actions may be taken to move in such directions.

- The problem addressed in this paper is then that of finding optimal exposures on a single underlying stock.

- Our approach differs from the existing literature in a several ways.
  - The focus here is on the instantaneous jump $x$ in the stock, which leads us to consider a possibly infinite Levy density rather than a probability measure
  - We assuming that the options that give rise to the exposure have a fixed maturity of 15 days, and that they will be unwound in 10 days. Thus, the base case scenario for the density of the instantaenous jump $x$ in the log price is the one calibrated from option prices with 5-day expiration.
  - We then maximize the expected exposure over all scenarios obtained by distorting the base density according to the MINMAXVAR distortion function $\Psi$.
  - As some of such scenarios are less likely than others, we add to the expected exposure under each scenarios a rebate that is zero for the base scenario and is higher for scenarios that are more far away from the base one.

- These considerations result in the following formulation of the maximization problem:

$$\begin{align}&\max_y \min_z\int_{\mathbb{R}}y(x)z(x)x^4m(x)dx-\frac{\int_{\mathbb{R}}y(x)m(x)x^2dx}{\int_{\mathbb{R}}m(x)x^2dx}\int_{\mathbb{R}}m(x)x^2dx+\int_{\mathbb{R}}|z(x)|^{\alpha}dx\\
    &\text{s.t. }z(x)\geq -\frac{1}{x^2}\\
    & \int_{\mathbb{R}}\left(z(x)x^2-(\lambda-1)\right)^+m(x)dx\leq \Phi(\lambda), \ \lambda\geq 1\\
    & \int_{\mathbb{R}}\left(-z(x)x^2-(1-\lambda)\right)^+m(x)dx\leq \tilde{\Phi}(\lambda), \ \lambda\leq 1\\
    & \int_{\mathbb{R}}y(x)\tilde{m}x^2dx=0\end{align}$$
  
- The densities $m$ and $\tilde{m}$ are the bilateral CGMY densities calibrated to the 5 and 15 day maturity options respectively.  

- The distortion and rebate parameters are set to specific values but may be obtained via backtesting.

- For a full description of the problem set up we refer to our accompanying paper.

## Formulation for discipline saddle programming

- To use the dsp extension of CVXPY we need to discretize the maximization problem

- Here we consider the problem

$$\begin{align}&\max_{y\in\mathbb{R}^K}\min_{z\in\mathbb{R}^N}z^TP_4\left(My-p_0^Ty\right)+p_2^T y +p_{2\alpha}^T |z|^{\alpha}\\
    & \text{ s.t. } z\geq 0\\
    & \qquad p_0^T\max(x_2^Tz-(\lambda-1),0)\leq \Phi(\lambda), \ \lambda\geq 1\\
    & \qquad p_0^T\max(1-\lambda-x_2^Tz,0)\leq -\tilde{\Phi}(\lambda), \ 0\leq \lambda\leq 1\\
    & \qquad q_0^TMy = 0.\end{align}$$

- where:
  - $P_k$ denotes the $N\times N$ diagonal matrix with diagonal elements given by $\delta m(x_j)x_j^2$, $j=1,...,2^k$,
  - $\mathbf{p}_{2\alpha}$, $\mathbf{p}_{2}$ are the Levy densities multiplied by $x^4$, $x^{2\alpha}$ and $x^2$ at the 5-day maturity
  - $\mathbf{p}_0$ and $\mathbf{q}_0$ are the Levy densities at the 5-day and 15-daymaturity
  - $\mathbf{y}$, $\mathbf{x}$ and $\mathbf{z}$ are vectors in $\mathbb{R}^{N}$
- This formulation of the problem satisfies the requirements to be solved via disciplined saddle programming, for which we refer to the github repository available at <https://github.com/cvxgrp/dsp> and the accompnying papers <https://arxiv.org/abs/2301.13427>
  - <https://arxiv.org/abs/2102.01002>.

## Additional Remarks

- The MATLAB folder contains a matlab m file, which solves the problem considered by running a python script

- As this only supports a previous version of Python, do not forget to match your vscode python version with that of the conda environment (v3.10.14)
