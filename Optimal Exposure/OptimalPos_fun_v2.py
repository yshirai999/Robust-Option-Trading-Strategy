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
    x: np.ndarray,
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

    
    z = cp.Variable(N)

#     y = cp.Variable(K)
#     f = dsp.inner( z, P @ ( MM @ y - q0 @ MM @ y) )
#     f1 = (p2) @ (MM @ y - q0 @ MM @ y)
    
    y = cp.Variable(1)
    f = dsp.inner( z, P @ ( cp.multiply(y,cp.exp(x)-1) - q0 @ cp.multiply(y,cp.exp(x)-1) ) )
    f1 = p2 @ ( cp.multiply(y,cp.exp(x)-1) - q0 @ cp.multiply(y,cp.exp(x)-1) )

    rho = pa2 @ cp.power(cp.abs(z),alpha)

    constraints = [z >= -x2inv]
    
    constraints.append(y>=-5000)
    constraints.append(y<=100)

    for i in range(Cl): 
        constraints.append(p0 @ cp.maximum( (1-lamn[i])-cp.multiply(x2,z), 0 ) <= -Phi_l[i])
    for i in range(Cu): 
        constraints.append(p0 @ cp.maximum( cp.multiply(x2,z)-(lamp[i]-1), 0 ) <= Phi_u[i])

    obj = dsp.MinimizeMaximize(rho+f+f1)
    prob = dsp.SaddlePointProblem(obj, constraints)
    prob.solve(solver = cp.CLARABEL, max_iter = 500, tol_infeas_abs = 1e-5, tol_feas = 1e-5, min_terminate_step_length = 1e-5, verbose = verbose)  # solves the problem
    #prob.solve(solver = cp.SCS, verbose = verbose)

    print(prob.value)
    return y.value

z = OptimalPos(pa2,p2,p4,p0,q0,MM,lamp,lamn,Phi_u,Phi_l,x,x2,x2inv,N,K,Cu,Cl,alpha,verbose)
