import dsp
import cvxpy as cp
import numpy as np
from Models import BG
from Distortions import minmaxvar as mmv
import matplotlib.pyplot as plt

p = np.array([2/3,1/6,1/6])
q = [1/3,1/3,1/3]

N = 100
a = np.linspace(0.5,2,N)

lam = 0.25
dist = mmv(lam)
Phi = dist.Phi(a)

theta = 0.75
alpha = 1.25
beta = 0.25
P = np.diag(p)
I = np.diag([1,1,1])

x = cp.Variable(3)
y = cp.Variable(3)
z = cp.Variable(3)
zx = cp.Variable(3)
zy = cp.Variable(3)

# # # print(np.shape(cp.vstack([P,P,P])))
# # #print(np.shape(cp.hstack([(x+y),-x,-y]).T))
# # #print(np.shape(cp.vstack([(x+y),-x,-y]).T))
# # print(np.shape(cp.vstack([z,zx,zy])))
# print(np.shape(cp.hstack([z,zx,zy])))

# # print(np.shape(P @ cp.vstack([(x+y),-x,-y]).T))

# print(np.shape(cp.reshape(P @ cp.vstack([(x+y),-x,-y]).T,(9,))))

# print(type((dsp.inner(cp.hstack([z,zx,zy]), np.reshape(P @ cp.vstack([(x+y),-x,-y]).T,(9,))))))

# print(np.shape((cp.multiply(theta, cp.power(z,alpha))+cp.multiply(1-theta,cp.power(z,-beta)))))
# f = (dsp.inner(cp.hstack([z,zx,zy]), np.reshape(P @ cp.vstack([(x+y),-x,-y]).T,(9,))))

f = dsp.inner(y, P @ z)

rho = p @ cp.multiply(theta, cp.power(z.T,alpha))+cp.multiply(1-theta,cp.power(z.T,-beta))

print(np.shape(cp.multiply(theta, cp.power(z.T,alpha))+cp.multiply(1-theta,cp.power(z.T,-beta))))
print(np.shape(rho))
print(np.shape(f))
# rhox = p @ (cp.multiply(theta, cp.power(zx,alpha))+cp.multiply(1-theta,cp.power(zx,-beta)))
# rhoy = p @ (cp.multiply(theta, cp.power(zy,alpha))+cp.multiply(1-theta,cp.power(zy,-beta)))
# obj = dsp.MinimizeMaximize(rho-rhox-rhoy-f)
obj = dsp.MinimizeMaximize(f+rho)
constraints = [p @ z == 1, z >= 0]
# for i in range(N): #range(len(a)):
#     constraints.append(p @ cp.maximum(z-a[i],0) <= Phi[i])
#     constraints.append(p @ cp.maximum(zx-a[i],0) <= Phi[i])
#     constraints.append(p @ cp.maximum(zy-a[i],0) <= Phi[i])

prob = dsp.SaddlePointProblem(obj, constraints)
prob.solve()  # solves the problem

print(prob.value)
