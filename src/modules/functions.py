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
    if(not param.cluster_mode):
        import matplotlib.pyplot as plt
        plt.figure()
        plt.plot(param.x,vx,label='The local potential')
        plt.grid()
        plt.legend()
        plt.show()
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
def occbkubkG_dns(param,occbk,ubkG):
    dns = np.zeros(param.NG,dtype='float64')
    work = np.empty_like(ubkG[:,0,0])
    NBact = np.shape(ubkG)[1]
    for ik in range(param.NK):
        for ib in range(NBact):
            work = np.fft.ifft(ubkG[:,ib,ik])
            dns = dns + occbk[ib,ik]*(np.abs(work))**2
    return dns

def occbkubkG_J(param,occbk,ubkG,A): #Exact formula should be checked=========
    J = 0.0
    for ik in range(param.NK):
        kpA = param.k[ik] + A
        for ib in range(param.NG):
            J = J + occbk[ib,ik]*(np.sum(param.G[:]*(np.abs(ubkG[:,ib,ik]))**2)*param.a/float(param.NG**2) + kpA)
    return J/param.a

def occbkubkG_Etot(param,occbk,ubkG,A): #Exact formula should be checked=========
    Etot = 0.0
    vx, vG, vGG, vGGk = get_vxvGvGGvGGk(param)
    hk = get_tGGk(param,0.0) + vGGk 
    for ik in range(param.NK):
        hubG = np.dot(hk[:,:,ik], ubkG[:,:,ik])
        for ib in range(param.NG):
            Etot = Etot + occbk[ib,ik]*np.real(np.vdot(ubkG[:,ib,ik],hubG[:,ib]))
    return Etot*param.a/float(param.NG**2) #orbital function is normalized to give correct number of particle in the cell.
#
def h_U(param,h):
    w, v = np.linalg.eigh(h)
    U = np.exp(-zI*w[0]*param.dt)*np.outer(v[0,:],np.conj(v[0,:])) + np.exp(-zI*w[1]*param.dt)*np.outer(v[1,:],np.conj(v[1,:]))
    return U
#
def Make_Efield(param):
    t = np.zeros([param.Nt],dtype=np.float64)
#    E = np.zeros([param.Nt],dtype=np.float64)
    E = np.zeros([param.Nt,param.Ncolor],dtype=np.float64)
    for it in range(param.Nt):
        t[it] = param.dt*it
    if (param.Ncolor == 1):
        icolor = 0
        for it in range(param.Nt):
            if (t[it] < param.Tpulse):
                E[it,icolor] = param.E0*(np.sin(pi*t[it]/param.Tpulse))**param.nenvelope*np.sin(param.omegac*(t[it] - 0.5*param.Tpulse) + param.phi_CEP)
    elif (param.Ncolor > 1):
        for icolor in range(param.Ncolor):
            for it in range(param.Nt):
                if (t[it] < param.Tpulse[icolor]):
                    E[it,icolor] = param.E0[icolor]*(np.sin(pi*t[it]/param.Tpulse[icolor]))**param.nenvelope[icolor]*np.sin(param.omegac[icolor]*(t[it] - 0.5*param.Tpulse[icolor]) + param.phi_CEP[icolor])
        E = np.sum(E,axis=1)
    else :
        print('ERROR: The parameter '+str(param.Ncolor)+' is improper.')
        sys.exit()
    return t, E
