import numpy as np

class minmaxvar():
    
    def __init__(self,
                lam: float
    ) -> None:
        self.lam = lam
    
    def Psi(self,
            x = np.ndarray
    ) -> np.ndarray:
        psi = 1 - ( 1-x**(1/(self.lam+1)) ) ** (1+self.lam)
        return psi

    def Phi(self,
            a = np.ndarray
    ) -> np.ndarray:
        phi = ( 1 / a ** (1/self.lam) ) ** (self.lam+1)
        return phi