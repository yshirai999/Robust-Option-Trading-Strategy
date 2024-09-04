import dsp
import cvxpy as cp
import numpy as np
import matplotlib.pyplot as plt

def OptimalPos(
    p_pa2: np.ndarray,
    p_na2: np.ndarray,
    p_p2: np.ndarray,
    p_n2: np.ndarray,
    p_p4: np.ndarray,
    p_n4: np.ndarray,
    A_pa2: np.ndarray,
    A_na2: np.ndarray,
    A_p2: np.ndarray,
    A_n2: np.ndarray,
    A_p4: np.ndarray,
    A_n4: np.ndarray,
    B_p: np.ndarray,
    p_p0: np.ndarray,
    B_n: np.ndarray,
    p_n0: np.ndarray,
    M: np.ndarray,
    lamp: np.ndarray,
    lamn: np.ndarray,
    Phi_u: np.ndarray,
    Phi_l: np.ndarray,
    x2: np.ndarray,
    x2inv: np.ndarray,
    N: float,
    K: float,
    C: float,
    theta: float = 0.25,
    alpha: float = 1.2,
    ):
    
    N = int(N)
    K = int(K)
    C = int(C)

    P = np.diag(A_p4*p_p4 + A_n4*p_n4)
    
    M = np.transpose(M)

    y = cp.Variable(K)
    zp = cp.Variable(N)
    zn = cp.Variable(N)
    
    f = dsp.inner( zp-zn, P @ ( M @ y ) )
    f1 = (A_p2*p_p2+A_n2*p_n2) @ (M @ y)
    rho = A_pa2*p_pa2 @ cp.multiply(theta, cp.power(zp+zn,alpha)) \
            + A_na2*p_na2 @ cp.multiply(theta, cp.power(zp+zn,alpha))

    constraints = [zp >= 0]
    constraints.append(zn >= 0)
    constraints.append(zn <= x2inv)
    for i in range(C): #range(len(a)):
        pp = B_p*p_p0
        pn = B_n*p_n0
        constraints.append(pp @ cp.maximum( cp.multiply(x2,cp.power(zp-zn,alpha))-(lamp[i]-1), 0 ) <= Phi_u[i])
        constraints.append(pn @ cp.maximum( (lamn[i]-1)-cp.multiply(x2,cp.power(zp-zn,alpha)), 0 ) <= -Phi_l[i])
#         constraints.append(pp @ cp.maximum( cp.multiply(x2,zp)-(lamp[i]-1), 0 ) <= Phi_u[i])
#         constraints.append(pn @ cp.maximum( (lamn[i]-1)-cp.multiply(x2,zp), 0 ) <= -Phi_l[i])

    obj = dsp.MinimizeMaximize(rho+f+f1)
    prob = dsp.SaddlePointProblem(obj, constraints)
    prob.solve(verbose = True)  # solves the problem

    print(prob.value)
    return y.value

z = OptimalPos(p_pa2,p_na2,p_p2,p_n2,p_p4,p_n4,
               A_pa2,A_na2,A_p2,A_n2,A_p4,A_n4,
               B_p,p_p0,B_n,p_n0,M,lamp,lamn,
               Phi_u,Phi_l,x2,x2inv,N,K,C,theta,alpha,)
