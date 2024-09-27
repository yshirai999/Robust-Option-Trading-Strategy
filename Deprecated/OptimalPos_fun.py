import dsp
import cvxpy as cp
import numpy as np
import matplotlib.pyplot as plt

def OptimalPos(
    pa2: np.ndarray,
    p2: np.ndarray,
    p4: np.ndarray,
    p0: np.ndarray,
    q0: np.ndarray,
    MM: np.ndarray,
    lamp: np.ndarray,
    lamn: np.ndarray,
    Phi_u: np.ndarray,
    Phi_l: np.ndarray,
    x2: np.ndarray,
    x2inv: np.ndarray,
    N: float,
    K: float,
    Cu: float,
    Cl: float,
    alpha: float,
    verbose: str
    ):
    
    N = int(N)
    K = int(K)
    Cu = int(Cu)
    Cl = int(Cl)
    
    verbose = verbose == 'True'

    P = np.diag(p4)

    MM = np.transpose(MM)

    y = cp.Variable(K)
    zp = cp.Variable(N)
    zn = cp.Variable(N)

    f = dsp.inner( zp-zn, P @ ( MM @ y - q0 @ MM @ y) )
    f1 = (p2) @ (MM @ y - q0 @ MM @ y)
    rho = pa2 @ cp.power(zp+zn,alpha)

    constraints = [zp >= 0]
    constraints.append(zn >= 0)
    constraints.append(zn <= x2inv)
    constraints.append(y>=-100)
    constraints.append(y<=100)

    for i in range(Cl): 
        constraints.append(p0 @ cp.maximum( (1-lamn[i])-cp.multiply(x2,zp-zn), 0 ) <= -Phi_l[i])
    for i in range(Cu): 
        constraints.append(p0 @ cp.maximum( cp.multiply(x2,zp-zn)-(lamp[i]-1), 0 ) <= Phi_u[i])

    obj = dsp.MinimizeMaximize(rho+f+f1)
    prob = dsp.SaddlePointProblem(obj, constraints)
    prob.solve(verbose = verbose)  # solves the problem

    print(prob.value)
    return y.value

z = OptimalPos(pa2,p2,p4,p0,q0,MM,lamp,lamn,Phi_u,Phi_l,x2,x2inv,N,K,Cu,Cl,alpha,verbose)
