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