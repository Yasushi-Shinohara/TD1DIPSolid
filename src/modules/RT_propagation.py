# coding: UTF-8
# This is created 2020/11/17 by Y. Shinohara
# This is lastly modified 2020/11/17 by Y. Shinohara
import sys
from modules.constants import tpi, zI
import numpy as np

class RT_propagation_class():
    def __init__(self):
        self.something = None

    def uGbk_forward(self, propagator_option):
        if (propagator_option.lower() == 'exp'):
            uGbk_forward = self.uGbk_forward_exp
            print('# The exponential expression for the temporal propagator is chosen.')
        elif (propagator_option.upper() == 'RK4'):
            uGbk_forward = self.uGbk_forward_RK4
            print('# The Runge-Kutta 4th for the temporal propagator is chosen.')
        elif ((propagator_option.upper() == 'RK4FFT') or (propagator_option.uppwer() == 'RK4_FFT')):
            uGbk_forward = self.uGbk_forward_RK4FFT
            print('# The Runge-Kutta 4th with FFT for the temporal propagator is chosen.')
        else :
            print('# ERROR: undefined propagator_option is called.')
            sys.exit()
        return uGbk_forward

    def hdt2U(self, hdt):
        eigs, coef = np.linalg.eigh(hdt)
        U = np.exp(-zI*eigs[0])*np.outer(coef[:,0],np.conj(coef[:,0]))
        for ib in range(1,len(eigs)):
            U = U + np.exp(-zI*eigs[ib])*np.outer(coef[:,ib],np.conj(coef[:,ib]))
        return U
        
    def uGbk_forward_exp(self, param, uGbk, hGGk, tGGk, vx):
        for ik in range(param.Nk):
            U = self.hdt2U(hGGk[:,:,ik]*param.dt)
            uGbk[:,:,ik] = np.dot(U, uGbk[:,:,ik])
        return uGbk

    def uh2hu(self, uGb, hGG):
        return np.dot(hGG,uGb)

    def uGbk_forward_RK4(self, param, uGbk, hGGk, tGGk, vx):
        for ik in range(param.Nk):
            k1 = self.uh2hu(uGbk[:,:,ik], hGGk[:,:,ik])/zI
            k2 = self.uh2hu(uGbk[:,:,ik] + 0.5*param.dt*k1, hGGk[:,:,ik])/zI
            k3 = self.uh2hu(uGbk[:,:,ik] + 0.5*param.dt*k2, hGGk[:,:,ik])/zI
            k4 = self.uh2hu(uGbk[:,:,ik] + param.dt*k3, hGGk[:,:,ik])/zI
            uGbk[:,:,ik] = uGbk[:,:,ik] + (k1 + 2.0*k2 + 2.0*k3 + k4)*param.dt/6.0 
        return uGbk

    def utv2hu_FFT(self, uGb, tGGdiag, vx):
        NBact = np.shape(uGb)[1]
        vuxb = np.empty_like(uGb)
        tuGb = np.empty_like(uGb)
        uxb = np.fft.ifft(uGb, axis=0)
        for ib in range(NBact):
            vuxb[:,ib] = vx[:]*uxb[:,ib]
            tuGb[:,ib] = tGGdiag[:]*uGb[:,ib]
        vuGb = np.fft.fft(vuxb, axis=0)
        return tuGb+vuGb

    def uGbk_forward_RK4FFT(self, param, uGbk, hGGk, tGGk, vx):
        tGGdiagk = np.diagonal(tGGk, axis1 = 0, axis2 = 1).T
        for ik in range(param.Nk):
            k1 = self.utv2hu_FFT(uGbk[:,:,ik], tGGdiagk[:,ik], vx)/zI
            k2 = self.utv2hu_FFT(uGbk[:,:,ik] + 0.5*param.dt*k1, tGGdiagk[:,ik], vx)/zI
            k3 = self.utv2hu_FFT(uGbk[:,:,ik] + 0.5*param.dt*k2, tGGdiagk[:,ik], vx)/zI
            k4 = self.utv2hu_FFT(uGbk[:,:,ik] + param.dt*k3, tGGdiagk[:,ik], vx)/zI
            uGbk[:,:,ik] = uGbk[:,:,ik] + (k1 + 2.0*k2 + 2.0*k3 + k4)*param.dt/6.0 
        return uGbk
