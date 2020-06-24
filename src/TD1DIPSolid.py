#!/usr/bin/python
# coding: UTF-8
# This is created 2020/06/25 by Y. Shinohara
# This is lastly modified 2020/06/25 by Y. Shinohara #This part is highly doubtable because of my lazyness
import time
ts = time.time()
from modules.print_funcs import print_header, print_footer, print_midtime, print_endtime
print_header()
import sys
import numpy as np
import math
import ctypes as ct
from modules.constants import *
from modules.parameters import parameter_class
param = parameter_class()
param.read_parameters()   #Initialization of the parameters and the replacement from the standard input
from modules.functions import get_hD, E_hOD, psih_Ene, h_U, Make_Efield
from modules.plot_funcs import plot_E, plot_RT

if (not param.cluster_mode):
    #Matplotlib is activated for the cluster_mode == True
    import matplotlib.pyplot as plt
    from matplotlib import cm #To include color map

#############################Prep. for the system########################
hD = get_hD(param)
hOD = E_hOD(param,0.0)
h = hD + hOD

psi = np.zeros([2],dtype=np.complex128)
psi[1] = 1.0  #Lower level is initailly fully occupied

Ene = psih_Ene(psi,h)
print('# System energy at initial:',Ene, '[a.u.] =',Ene*Hartree, ' [eV]')
#sys.exit()

#############################Prep. for RT################################
t, E = Make_Efield(param)
if (param.PC_option):
    Eave = 0.0*E
    for it in range(param.Nt-1):
        Eave[it] = 0.5*(E[it] + E[it+1])
    Eave[param.Nt - 1] = 1.0*Eave[param.Nt - 2]
nv = np.zeros([param.Nt],dtype=np.float64)
nc = np.zeros([param.Nt],dtype=np.float64)
Ene = np.zeros([param.Nt],dtype=np.float64)

if (np.amax(t) < np.amax(param.Tpulse)):
    print('# WARNING: max(t) is shorter than Tpulse')
        
if (not param.cluster_mode):
    #Plot shape of the electric field
    plot_E(plt,cm, t,E)

tt = time.time()
print_midtime(ts,tt)
#############################RT calculation##############################
#Time-propagation
for it in range(param.Nt):
    if (param.PC_option):
        hOD = E_hOD(param,Eave[it])
    else:
        hOD = E_hOD(param,E[it])
    h = hD + hOD
    U = h_U(param,h)
    psi = np.dot(U, psi)
    nv[it] = (np.abs(psi[0]))**2
    nc[it] = (np.abs(psi[1]))**2
    norm = np.linalg.norm(psi)
    Ene[it] = psih_Ene(psi,h)
    if (it%1000 == 0):
        print('# ',it, Ene[it], norm)

print('# System energy at end:',Ene[param.Nt-1], '[a.u.] =',Ene[param.Nt-1]*Hartree, ' [eV]')
print('# Absorbed energy:',Ene[param.Nt-1]-Ene[0], '[a.u.] =',(Ene[param.Nt-1]-Ene[0])*Hartree, ' [eV]')

if (not param.cluster_mode):
    #Plot data obtained in real-time evolution, nv, nc, Ene
    plot_RT(plt,cm, t,nv,nc,Ene)

te = time.time()
print_endtime(ts,tt,te,param.Nt)

if (param.minimal_output):
    np.savez('tEne.npz', t=t, Ene=Ene)
else:
    np.savez('RTall.npz', t=t, nv=nv, nc=nc, Ene=Ene)

print_footer() 
sys.exit()

