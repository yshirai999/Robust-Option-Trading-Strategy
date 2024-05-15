import dsp
import cvxpy as cp
import numpy as np
from Models import BG
from Distortions import minmaxvar as mmv
import matplotlib.pyplot as plt

params = [0.04,10/52,1,5/52,0.04,10/52,1,5/52]

k = 8
W = 1
M = [-1,1]
model = BG(k,M,params)
p = model.p
q = model.q
x = model.y

# fig = plt.figure()
# axes = fig.add_axes([0.1, 0.1, 1, 1])
# axes.set_xlim(M[0], M[1])
# axes.set_ylim(0, max([max(p),max(q)]))
# axes.plot(x, p)
# axes.plot(x, q)
# plt.show()

lam = 5
a = np.linspace(0.1,5,100)
dist = mmv(lam)
Phi = dist.Phi(a)
# fig = plt.figure()
# axes = fig.add_axes([0.1, 0.1, 1, 1])
# axes.set_xlim(a[0]-0.5, a[-1]+0.5)
# axes.set_ylim(0, max(Phi))
# axes.plot(a, Phi)
# # axes.plot(a,dist.Psi(a))
# plt.show()

theta = 0.75
alpha = 1.25
beta = 0.25
P = np.diag(p)

y = cp.Variable(2**k)
z = cp.Variable(2**k)

f = dsp.inner(z, P @ (y - cp.multiply(W,cp.exp(x))))
rho = p @ (cp.multiply(theta, cp.power(z,alpha))+cp.multiply(1-theta,cp.power(z,-beta)))
obj = dsp.MinimizeMaximize(rho-f)
constraints = [q @ y == W, p @ z == 1, z >= 0]
for i in range(89): #range(len(a)):
    constraints.append(p @ cp.maximum(z-a[i],0) <= Phi[i])

'''
Note: too many constraints generate some problems in CVXPY, see https://github.com/cvxpy/cvxpy/issues/826
Because of this I have here limited the number of constraints. It can be seen that the solution z does satisfy the constraint.

A potential solution is suggested on github:
The problem is that when calculating the shape for scipy.sparse.csc_matrix, it always gets a negative value,
since the shape is numpy.int32 and if the shape is too large, we get overflow problem. The idea to solve it is very simple,
just convert it into numpy.int64.
So do the following steps:
1. find canonInterface.py in cvxpy\cvxcore\python and find get_problem_matrix function. you will see around 368line, use this
instead A = scipy.sparse.csc_matrix((V, (I, J)), shape=(np.int64(constr_length)*np.int64(var_length+1), param_size_plus_one))
2. find conic_solver.py in \cvxpy\reductions\solvers\conic_solvers and find format_constraints function. 
You will see around 225 line and use this instead:
restructured_A = 
    restructured_A.reshape(np.int64(restruct_mat.shape[0]) * np.int64(problem.x.size + 1),problem.A.shape[1], order='F')
'''

obj = dsp.MinimizeMaximize(f)
prob = dsp.SaddlePointProblem(obj, constraints)
prob.solve()  # solves the problem

print(prob.value)

fig = plt.figure()
axes = fig.add_axes([0.1, 0.1, 0.75, 0.75])
axes.set_xlim(np.log(W)+M[0], np.log(W) + M[-1])
axes.set_ylim(min(y.value), max(y.value))
axes.plot(x,y.value)
# axes.plot(a,dist.Psi(a))
plt.show()

# Constraint satisfied
# print(q @ y.value)
zz = z.value
N = len(a)
const = []
for i in range(N): #range(len(a)):
    const.append(p @ np.maximum(zz-a[i],0))

const = np.array(const)
fig = plt.figure()
axes = fig.add_axes([0.1, 0.1, 0.75, 0.75])
axes.set_xlim(a[0], a[-1])
axes.set_ylim(min(const), max(Phi))
axes.plot(a,const)
axes.plot(a,Phi)
# axes.plot(a,dist.Psi(a))
plt.show()