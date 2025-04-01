import dsp
import cvxpy as cp
import numpy as np
import matplotlib.pyplot as plt

def OptimalPos(
    p: np.ndarray,
    q: np.ndarray,
    x: np.ndarray,
    a: np.ndarray,
    Phi: np.ndarray,
    k: float = 8,
    theta: float = 0.25,
    alpha: float = 0.25,
    beta: float = 0.25,
    W: float = 1,
    N: float = 50
    ):
    
    k = int(k)
    N = int(N)
    
    P = np.diag(p)

    y = cp.Variable(2**k)
    z = cp.Variable(2**k)

    f = dsp.inner(z, P @ (y - q @ y))
    rho = p @ (cp.multiply(theta, cp.power(z,alpha))+cp.multiply(1-theta,cp.power(z,-beta)))
    obj = dsp.MinimizeMaximize(rho+f)
    constraints = [p @ z == 1, z >= 0]
    for i in range(N): #range(len(a)):
        constraints.append(p @ cp.maximum(z-a[i],0) <= Phi[i])

    prob = dsp.SaddlePointProblem(obj, constraints)
    prob.solve()  # solves the problem

    print(prob.value)

    # fig = plt.figure()
    # axes = fig.add_axes([0.1, 0.1, 0.75, 0.75])
    # M = [-0.3,0.3]
    # axes.set_xlim(np.log(W)+M[0], np.log(W)+M[-1])
    # #axes.set_ylim(min(y.value), max(y.value))
    # #axes.plot(x,y.value-W*np.exp(x)
    # axes.plot(x,y.value)
    
    # plt.savefig('Pos.png')
    # plt.show()
        

    # # Constraint satisfied
    # print(q @ y.value)
    # zz = z.value
    # N = len(a)
    # const = []
    # for i in range(N): #range(len(a)):
    #     const.append(p @ np.maximum(zz-a[i],0))

    # const = np.array(const)
    # fig = plt.figure()
    # axes = fig.add_axes([0.1, 0.1, 0.75, 0.75])
    # axes.set_xlim(a[0], a[-1])
    # axes.set_ylim(min(const), max(Phi))
    # axes.plot(a,const)
    # axes.plot(a,Phi)
    # # axes.plot(a,dist.Psi(a))
    # plt.show()

    return y.value

z = OptimalPos(p, q, x, a, Phi, k, theta, alpha, beta, W, N)
