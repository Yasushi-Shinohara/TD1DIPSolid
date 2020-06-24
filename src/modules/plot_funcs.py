#!/usr/bin/python
# coding: UTF-8
# This is created 2020/04/20 by Y. Shinohara
# This is lastly modified 2020/04/20 by Y. Shinohara
from modules.constants import *

def plot_E(plt,cm, t,E):
    plt.figure()
    plt.title('Electric field')
    plt.xlabel('Time [fs]')
    plt.ylabel('Field strength [V/nm]')
    plt.xlim(0.0,np.amax(t)*Atomtime)
    plt.plot(t*Atomtime,E*Atomfield)
    plt.grid()
    plt.show()
#
def plot_RT(plt,cm, t,nv,nc,Ene):
    plt.figure()
    plt.title('Population of the lower and the upper level')
    plt.xlabel('Time [fs]')
    plt.ylabel('Populations')
    plt.xlim(0.0,np.amax(t)*Atomtime)
    plt.plot(t*Atomtime,nv,label='nv')
    plt.plot(t*Atomtime,nc,label='nc')
    plt.legend()
    plt.grid()
    plt.show()
#
    plt.figure()
    plt.title('Population of the lower level with log-scale')
    plt.yscale('log')
    plt.xlabel('Time [fs]')
    plt.ylabel('Populations')
    plt.xlim(0.0,np.amax(t)*Atomtime)
    plt.plot(t*Atomtime,nv,label='nv')
    plt.legend()
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
