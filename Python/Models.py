import numpy as np
import scipy as sp
import scipy.special as sc

class BG():

    def __init__(self,
                k: int,
                M: list[float],
                params: list[float]
    ):
        self.k = k
        self.M = M
        self.y = np.linspace(M[0],M[1],2**k)
        ip = (self.y > 0)
        im = (self.y < 0)
        yp = self.y[ip]
        ym = -self.y[im]
        self.delta = (M[1]-M[0])/(2**k-1)

        self.bp = params[0]
        self.cp = params[1]
        self.bn = params[2]
        self.cn = params[3]            
        b = (1/self.bp+1/self.bn)
        gammap = sc.gamma(self.cp)
        gamman = sc.gamma(self.cn)
        lam = 0.5*(self.cp+self.cn)
        mu = 0.5*(self.cp+self.cn-1)
        x = self.y*b
        xp = x[ip]
        xm = -x[im]
        
        W = np.zeros(x.shape)
        p = np.zeros(x.shape)
        W[ip] = np.exp(-0.5*xp) * (xp**(mu+0.5)) * sc.hyperu(mu-lam+0.5, 1+2*mu, xp)
        W[im] = np.exp(-0.5*xm) * (xm**(mu+0.5)) * sc.hyperu(mu-lam+0.5, 1+2*mu, xm)
        p[ip] = ( (self.bp)**(-self.cp) ) * ( (self.bn)**(-self.cn) ) * ( (yp)**(0.5*(self.cp+self.cn)-1) ) * np.exp(-0.5*yp) * W[ip] / gammap
        p[im] = ( (self.bp)**(-self.cp) ) * ( (self.bn)**(-self.cn) ) * ( (ym)**(0.5*(self.cp+self.cn)-1) ) * np.exp(-0.5*ym) * W[im] / gamman
        self.p = p * self.delta / ( b**(0.5*(self.cp+self.cn)) )
        self.p = p/sum(p)

        self.bptil = params[4]
        self.cptil = params[5]
        self.bntil = params[6]
        self.cntil = params[7]            
        b = (1/self.bptil+1/self.bntil)
        gammap = sc.gamma(self.cptil)
        gamman = sc.gamma(self.cntil)
        lam = 0.5*(self.cptil+self.cntil)
        mu = 0.5*(self.cptil+self.cntil-1)
        x = self.y*b
        ip = (x > 0)
        im = (x < 0)
        xp = x[ip]
        xm = -x[im]
        W = np.zeros(x.shape)
        q = np.zeros(x.shape)
        W[ip] = np.exp(-0.5*xp) * (xp**(mu+0.5)) * sc.hyperu(mu-lam+0.5, 1+2*mu, xp)
        W[im] = np.exp(-0.5*xm) * (xm**(mu+0.5)) * sc.hyperu(mu-lam+0.5, 1+2*mu, xm)
        q[ip] = ( (self.bptil)**(-self.cptil) ) * ( (self.bntil)**(-self.cntil) ) * ( (yp)**(0.5*(self.cptil+self.cntil)-1) ) * np.exp(-0.5*yp) * W[ip] / gammap
        q[im] = ( (self.bptil)**(-self.cptil) ) * ( (self.bntil)**(-self.cntil) ) * ( (ym)**(0.5*(self.cptil+self.cntil)-1) ) * np.exp(-0.5*ym) * W[im] / gamman
        q = q * self.delta / ( b**(0.5*(self.cptil+self.cntil)) )
        self.q = q/sum(q)

    def BGprice(self,
                S: list[float],
                k0: np.ndarray,
                k1: np.ndarray,
                r: float,
                T: float,
                alpha: float, lam: float, eta: float, N: int,
    ) -> list[np.ndarray]:
        bp = self.bptil
        bn = self.bntil
        cp = self.cptil
        cn = self.bntil

        beta = np.log(S[0])-lam*N/2
        k = beta+(np.cumsum(np.ones(N,1))-1)*lam
        u = (np.cumsum(np.ones(N,1))-1)*eta
        w = np.ones(N,1)*eta
        w[0] = w[0]/2
        x = np.exp(-1j*beta*u)*self.Psi_BG(u,bp,cp,bn,cn,T,r,alpha,S[0])*w
        Call = np.real((np.exp(-alpha*k)/np.pi)*sp.fft(x,N))
        kk = np.log(k0)
        C0 = np.interp(k,Call,kk)
        P0 = C0 - S[0] + k0*np.exp(-r*T)

        beta = np.log(S[1])-lam*N/2
        k = beta+(np.cumsum(np.ones(N,1))-1)*lam
        u = (np.cumsum(np.ones(N,1))-1)*eta
        w = np.ones(N,1)*eta
        w[0] = w[0]/2
        x = np.exp(-1j*beta*u)*self.Psi_BG(u,bp,cp,bn,cn,T,r,alpha,S[0])*w
        Call = np.real((np.exp(-alpha*k)/np.pi)*sp.fft(x,N))
        kk = np.log(k1)
        C1 = np.interp(k,Call,kk)
        P1 = C1 - S[1] + k1*np.exp(-r*T)

        return [C0, P0, C1, P1]

    def Phi_BG(self,
                u: float,
                r: float,
                k: np.ndarray,
                T: float
    ) -> list[float]:
        bp = self.bptil
        bn = self.bntil
        cp = self.cptil
        cn = self.bntil
        d0 = 1j*u*(np.log(k)+T*(r+cp*np.log(1-bp)+cn*np.log(1+bn)))
        Phi = np.exp(d0)*((1-1j*u*bp)^(-T*cp))*((1+1j*u*bn)^(-T*cn))
        return Phi

    def Psi_BG(self,
               u: float,
               r: float,
               k: np.ndarray,
               T: float,
               alpha: float
        ) -> list[float]: 
        bp = self.bptil
        bn = self.bntil
        cp = self.cptil
        cn = self.bntil
        Psi = (np.exp(-r*T)/((alpha+1j*u)/(alpha+1j*u+1)))/self.Phi_BG(u-(alpha+1)*1j,bp,cp,bn,cn,T,k,r)
        return Psi