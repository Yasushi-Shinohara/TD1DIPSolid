# coding: UTF-8
# This is created 2020/11/17 by Y. Shinohara
# This is lastly modified 2020/11/17 by Y. Shinohara
import sys
from modules.constants import tpi, zI
import numpy as np

class RT_propagation_class():
    def __init__(self):
        self.something = None

    def hdt2U(self, hdt):
        eigs, coef = np.linalg.eigh(hdt)
        U = np.exp(-zI*eigs[0])*np.outer(coef[:,0],np.conj(coef[:,0]))
        for ib in range(1,len(eigs)):
            U = U + np.exp(-zI*eigs[ib])*np.outer(coef[:,ib],np.conj(coef[:,ib]))
        return U
        
    def uGbkhGGk2uGbk_exp(self, param, uGbk, hGGk):
        for ik in range(param.Nk):
            U = self.hdt2U(hGGk[:,:,ik]*param.dt)
            uGbk[:,:,ik] = np.dot(U, uGbk[:,:,ik])
        return uGbk

    def uh2hu(self, param, uGb, hGG):
        return np.dot(hGG,uGb)

    def uGbkhGGk2uGbk_RK4(self, param, uGbk, hGGk):
        for ik in range(param.Nk):
            k1 = self.uh2hu(self, uGbk[:,:,ik], hGGk[:,:,ik])/zI
            k2 = self.uh2hu(self, uGbk[:,:,ik] + 0.5*param.dt*k1, hGGk[:,:,ik])/zI
            k3 = self.uh2hu(self, uGbk[:,:,ik] + 0.5*param.dt*k2, hGGk[:,:,ik])/zI
            k4 = self.uh2hu(self, uGbk[:,:,ik] + param.dt*k3, hGGk[:,:,ik])/zI
            uGbk[:,:,ik] = uGbk[:,:,ik] + (k1 + 2.0*k2 + 2.0*k3 + k4)*param.dt/6.0 
        return uGbk
