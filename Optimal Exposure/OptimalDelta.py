import dsp
import cvxpy as cp
import numpy as np

def OptimalDelta(
    wwa: np.ndarray,
    ww1: np.ndarray,
    ww2: np.ndarray,
    ww0: np.ndarray,
    wwq: np.ndarray,
    lamp: np.ndarray,
    lamn: np.ndarray,
    Phip: np.ndarray,
    Phin: np.ndarray,
    A: np.ndarray,
    xxxxx: np.ndarray,
    xxxxx2: np.ndarray,
    xxxxx2b: np.ndarray,
    alpha: float,
    ky: float,
    kz: float,
    nlp: float,
    nln: float,
    verbose: str
    ):
    ky = int(ky)
    kz = int(kz)
    
   

    nlp = int(nlp)
    nln=int(nln)
    
    verbose = verbose == 'True'

    P = np.diag(ww2)
  #  Q = np.diag(ww1)
    A = np.transpose(A)
    y = cp.Variable(ky)
    z = cp.Variable(kz)
    

    costy = ww1 @ (A @ y)
    costnum=ww0 @ xxxxx2
    posnum=costy/costnum

    #g = posnum * (ww2 @ z)

    #g=dsp.inner(z, P @ Q @ (A @ y) )
  
    f1 = ww1 @ (A @ y)

    f = dsp.inner( z, P @ ( (A @ y) - ww0 @ (A @ y)/costnum) )
    rho = wwa @ cp.power(cp.abs(z),alpha)
    # rho=p @ (cp.multiply(1/(1+theta), cp.power(z,(1+theta)))+cp.multiply(1/(1-theta),cp.power(z,(1-theta)))) + 2*theta/(1-theta**2)
    constraints = [ z >= -xxxxx2b,  y >= -5000, y <= 5000, wwq @ (A @ y)==0]
    for i in range(nlp):
            constraints.append(ww0 @ cp.maximum(cp.multiply(z,xxxxx2) -(lamp[i]-1),0) <= Phip[i])
    for i in range(nln):
            constraints.append(ww0 @ cp.maximum(cp.multiply(-z,xxxxx2) - (1-lamn[i]),0) <= Phin[i]) 
    obj=dsp.MinimizeMaximize(rho+f+f1)                 
    prob=dsp.SaddlePointProblem(obj,constraints)
    prob.solve(solver = cp.CLARABEL, max_iter = 500, tol_infeas_abs = 1e-5, tol_feas = 1e-5, min_terminate_step_length = 1e-5, verbose = verbose)

    return y.value , z.value

q = OptimalDelta(wwa,ww1,ww2,ww0,wwq,lamp,lamn,Phip,Phin,A,xxxxx,xxxxx2,xxxxx2b,alpha,ky,kz,nlp,nln,verbose)