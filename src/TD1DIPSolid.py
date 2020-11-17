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
#DEBUGDEBUG
param.read_parameters()    #Initialization of the parameters and the replacement from the standard input#
#DEBUGDEBUG
param.grid_constructions() #
param.get_Nocc()           #
from modules.functions import * #Caution!! This should be modified 
from modules.plot_funcs import plot_E, plot_RT

if (not param.cluster_mode):
    #Matplotlib is activated for the cluster_mode == True
    import matplotlib.pyplot as plt
    from matplotlib import cm #To include color map

#############################Prep. for the system########################
uGbk = np.zeros([param.NG, param.NG, param.NK],dtype='complex128') #Wave function in reciprocal space
hGGk = np.zeros([param.NG, param.NG, param.NK],dtype='complex128') #Hamiltonian in terms of reciprocal space
epsbk = np.zeros([param.NG, param.NK],dtype='float64') #Hamiltonian in terms of reciprocal space
occbk = np.zeros([param.NG, param.NK],dtype='float64') #Occupation number
occbk[0:param.Nocc,:] = 2.0/float(param.NK)

vx, vG, vGG, vGGk = get_vxvGvGGvGGk(param)
tGGk = get_tGGk(param,0.0)
hGGk = tGGk + vGGk

#Band calculation 
for ik in range(param.NK):
    epsbk[:,ik], uGbk[:,:,ik] = np.linalg.eigh(hGGk[:,:,ik])
uGbk = uGbk/np.sqrt(param.a)*float(param.NG) #Normalization
Eg = np.amin(epsbk[param.Nocc,:])-np.amax(epsbk[param.Nocc - 1,:])
print('Eg = '+str(Eg)+' a.u. = '+str(Hartree*Eg)+' eV')

if (not param.cluster_mode):
    plt.figure()
    plt.xlim(np.amin(param.k), np.amax(param.k))
    plt.xlabel('$k$ [a.u.]')
    plt.ylabel('$\epsilon_{bk}$ [eV]')
    for ib in range(4):
        plt.plot(param.k,epsbk[ib,:]*Hartree)
    plt.grid()
    plt.show()

print('Band calculation is done properly.    ')
print('######################################')

#Ene = psih_Ene(psi,h)
#Ene = 
dns = occbkuGbk_dns(param,occbk,uGbk)
print('## Check for dns, '+str(np.sum(dns)*param.H))
J = occbkuGbk_J(param,occbk,uGbk,0.0) #Matter current
print('## Check for current, '+str(J))
Ene = occbkuGbkhGGk_Ene(param,occbk,uGbk,hGGk)
print('## Check for Ene, '+str(Ene))
print('# System energy at initial:',Ene, '[a.u.] =',Ene*Hartree, ' [eV]')

#sys.exit()

#############################Prep. for RT################################
t, A, E = Make_Efield(param)
if (param.PC_option):
    Eave = 0.0*E
    Aave = 0.0*A
    for it in range(param.Nt-1):
        Eave[it] = 0.5*(E[it] + E[it+1])
        Aave[it] = 0.5*(A[it] + A[it+1])
    Eave[param.Nt - 1] = 1.0*Eave[param.Nt - 2]
    Aave[param.Nt - 1] = 1.0*Aave[param.Nt - 2]
nv = np.zeros([param.Nt],dtype=np.float64)
nc = np.zeros([param.Nt],dtype=np.float64)
Ene = np.zeros([param.Nt],dtype=np.float64)
J = np.zeros([param.Nt],dtype=np.float64)

if (np.amax(t) < np.amax(param.Tpulse)):
    print('# WARNING: max(t) is shorter than Tpulse')
        
if (not param.cluster_mode):
    #Plot shape of the electric field
    plot_E(plt,cm, t,E)

tt = time.time()
print_midtime(ts,tt)
#sys.exit()
#############################RT calculation##############################
#Time-propagation
from modules.RT_propagation import RT_propagation_class
from modules.parameters import parameter_class
RTc = RT_propagation_class()
RT_option = 'exp'
if (RT_option == 'exp'):
    uGbkhkGG2uGbk = RTc.uGbkhGGk2uGbk_exp
else :
    print('# ERROR: undefined RT_option is called.')
    sys.exit()

for it in range(param.Nt):
    J[it] = occbkuGbk_J(param,occbk,uGbk,A[it])
    Ene[it] = occbkuGbkhGGk_Ene(param,occbk,uGbk,hGGk)
    if (param.PC_option):
        tGGk = get_tGGk(param,Aave[it])
    else:
        tGGk = get_tGGk(param,A[it])
    hGGk = tGGk + vGGk
    uGbk = uGbkhkGG2uGbk(param, uGbk, hGGk)
    if (it%1000 == 0):
        dns = occbkuGbk_dns(param,occbk,uGbk)
        print(it,np.sum(dns)*param.H, J[it], Ene[it])

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

