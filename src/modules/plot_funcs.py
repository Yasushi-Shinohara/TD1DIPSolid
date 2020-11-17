#!/usr/bin/python
# coding: UTF-8
# This is created 2020/04/20 by Y. Shinohara
# This is lastly modified 2020/04/20 by Y. Shinohara
from modules.constants import *

def plot_band(plt,cm, param, epsbk, Nbandmax = 4):
    plt.figure()
    plt.xlim(np.amin(param.k), np.amax(param.k))
    plt.xlabel('$k$ [a.u.]')
    plt.ylabel('$\epsilon_{bk}$ [eV]')
    for ib in range(Nbandmax):
        plt.plot(param.k,epsbk[ib,:]*Hartree)
    plt.grid()
    plt.show()
#
def plot_AE(plt,cm, param, t, A, E):
    plt.figure()
    plt.title('Vector potential')
    plt.xlabel('Time [fs]')
    plt.ylabel('Field strength [V/nm]')
    plt.xlim(0.0,np.amax(t)*Atomtime)
    plt.plot(t*Atomtime,A)
    plt.plot(t*Atomtime,param.b*np.ones(param.Nt)/2.0)
    plt.plot(t*Atomtime,-param.b*np.ones(param.Nt)/2.0)
    plt.grid()
    plt.show()
#
    plt.figure()
    plt.title('Electric field')
    plt.xlabel('Time [fs]')
    plt.ylabel('Field strength [V/nm]')
    plt.xlim(0.0,np.amax(t)*Atomtime)
    plt.plot(t*Atomtime,E*Atomfield)
    plt.grid()
    plt.show()
#
def plot_RT(plt,cm, t,J,Ene):
    plt.figure()
    plt.xlabel('Time [fs]')
    plt.ylabel('Current density [a.u.]')
    plt.xlim(0.0,np.amax(t)*Atomtime)
    plt.plot(t*Atomtime,J)
    plt.grid()
    plt.show()
#
    plt.figure()
    plt.title('Energy of the system')
    plt.xlabel('Time [fs]')
    plt.ylabel('Energy [eV]')
    plt.xlim(0.0,np.amax(t)*Atomtime)
    plt.plot(t*Atomtime,Ene*Hartree)
    plt.grid()
    plt.show()
