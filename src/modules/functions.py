# coding: UTF-8
#Relevant functions are written
# This is created 2020/04/17 by Y. Shinohara
# This is lastly modified 2020/05/20 by Y. Shinohara
import os
import math
import numpy as np
from modules.constants import *
#
def get_vxvGvGGvGGk(param):
    alpha = 5.0e-2
    beta = 5.0e-2
    gamma = 1.0e-1
    v0 =  0.37
    vx = -v0*(1.0+np.cos(tpi*param.x/param.a))
    vG = np.fft.fft(vx)/np.float(param.NG)
    vGG = np.zeros([param.NG, param.NG], dtype='complex128') 
    for ig1 in range(param.NG):
        gind = np.remainder(ig1 - np.arange(param.NG), param.NG)
        for ig2 in range(param.NG):
            igloc = gind[ig2]
            vGG[ig1,ig2] = vG[igloc] #This definition should be carefully checked 
    vGGk = np.zeros([param.NG, param.NG, param.NK], dtype='complex128')
    for ik in range(param.NK):
        vGGk[:,:,ik] = 1.0*vGG[:,:]
    return vx, vG, vGG, vGGk
#
def get_tGGk(param, A):
    tGGk = np.zeros([param.NG, param.NG, param.NK],dtype='complex128')
    for ik in range(param.NK):
        kpA = param.k[ik] + A
        for ig in range(param.NG):
            tGGk[ig, ig ,ik] = tGGk[ig, ig, ik] + 0.5*(param.G[ig] + kpA)**2
    return tGGk

#Relevant functions
def occbkuGbk_dns(param,occbk,uGbk):
    dns = np.zeros(param.NG,dtype='float64')
    work = np.empty_like(uGbk[:,0,0])
    NBact = np.shape(uGbk)[1]
    for ik in range(param.NK):
        for ib in range(NBact):
            work = np.fft.ifft(uGbk[:,ib,ik])
            dns = dns + occbk[ib,ik]*(np.abs(work))**2
    return dns

def occbkuGbk_J(param,occbk,uGbk,A): #Exact formula should be checked=========
    J = 0.0
    for ik in range(param.NK):
        kpA = param.k[ik] + A
        for ib in range(param.NG):
            J = J + occbk[ib,ik]*(np.sum(param.G[:]*(np.abs(uGbk[:,ib,ik]))**2)*param.a/float(param.NG**2) + kpA)
    return J/param.a

def occbkuGbkhGGk_Ene(param,occbk,uGbk, hGGk): #Exact formula should be checked=========
    Ene = 0.0
    hk = 1.0*hGGk
    for ik in range(param.NK):
        hubG = np.dot(hk[:,:,ik], uGbk[:,:,ik])
        for ib in range(param.NG):
            Ene = Ene + occbk[ib,ik]*np.real(np.vdot(uGbk[:,ib,ik],hubG[:,ib]))
    return Ene*param.a/float(param.NG**2) #orbital function is normalized to give correct number of particle in the cell.
#

def Make_Efield(param):
    t = np.zeros([param.Nt],dtype=np.float64)
    E = np.zeros([param.Nt],dtype=np.float64)
    A = np.zeros([param.Nt,param.Ncolor],dtype=np.float64)
    for it in range(param.Nt):
        t[it] = param.dt*it
    if (param.Ncolor == 1):
        icolor = 0
        for it in range(param.Nt):
            if (t[it] < param.Tpulse):
                A[it,icolor] = (param.E0/param.omegac)*(np.sin(pi*t[it]/param.Tpulse))**param.nenvelope*np.cos(param.omegac*(t[it] - 0.5*param.Tpulse) + param.phi_CEP)
    elif (param.Ncolor > 1):
        for icolor in range(param.Ncolor):
            for it in range(param.Nt):
                if (t[it] < param.Tpulse[icolor]):
                    A[it,icolor] = (param.E0[icolor]/param.omegac[icolor])*(np.sin(pi*t[it]/param.Tpulse[icolor]))**param.nenvelope[icolor]*np.cos(param.omegac[icolor]*(t[it] - 0.5*param.Tpulse[icolor]) + param.phi_CEP[icolor])
        A = np.sum(E,axis=1)
    else :
        print('ERROR: The parameter '+str(param.Ncolor)+' is improper.')
        sys.exit()

    for it in range(1,param.Nt-1):
        E[it] = -(A[it+1] - A[it-1])/2.0/param.dt
    E[0] = 2.0*E[1] - E[2]
    E[param.Nt-1] = 2.0*E[param.Nt-2] - E[param.Nt-3]    
    return t, A, E
