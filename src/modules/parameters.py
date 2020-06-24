# coding: UTF-8
# This is created 2020/04/17 by Y. Shinohara
# This is lastly modified 2020/05/20 by Y. Shinohara
import sys
from modules.constants import tpi, Atomtime, Hartree, Atomfield, aB
import numpy as np

class parameter_class:
    def __init__(self):
        #Default values for the parameters
        self.cluster_mode = False
        self.PC_option = True         #Predictor-corrector option
        self.minimal_output = True    #A flag to write-out minimal data or not
        self.a = 10.0                 #A spatial dimension constant that the dipole : d = e*a
        self.Delta = 9.0/Hartree      #The Gap of the two-level system in a.u.
        self.dt = 5.0e-1              #The size of the time-step
        self.Nt = 4000                #The number of time steps
        self.Ncolor = 1               #Number of color for the field
        self.omegac = 1.55/Hartree    #Photon energy
        self.phi_CEP = 0.0/tpi        #Carrier envelope phase
        self.Tpulse = 40.0/Atomtime   #A parameter for pulse duration
        self.nenvelope = 4            #Power for the sing envelope function
        self.E0 = 1.0e0/Atomfield     #Field strength in a.u.

    def read_parameters(self):
        argv = sys.argv
        argc = len(argv)
        #Reading data from the standard input
        if (argc == 1):
            print('# The default parameters are chosen.')
        elif (argc == 2):
            print('# Name of input file is "'+argv[1]+'".')
            f = open(argv[1],'r')
            lines = f.readlines()
            f.close
            Nlen = len(lines)
            text = [0]*Nlen
            for i in range(Nlen):
                text[i] = lines[i].strip()
            for i in range(Nlen):
                if (str(text[i]) == 'cluster_mode'):
                    if (str(text[i+1]) == 'True'):
                        self.cluster_mode = True
                    else:
                        self.cluster_mode = False
                if (str(text[i]) == 'PC_option'):
                    if (str(text[i+1]) == 'True'):
                        self.PC_option = True
                    else:
                        self.PC_option = False
                if (str(text[i]) == 'minimal_output'):
                    if (str(text[i+1]) == 'True'):
                        self.minimal_output = True
                    else:
                        self.minimal_output = False
                if (str(text[i]) == 'a'):
                    self.a = float(str(text[i+1]))
                if (str(text[i]) == 'Delta'):
                    self.Delta = float(str(text[i+1]))
                if (str(text[i]) == 'dt'):
                    self.dt = float(str(text[i+1]))
                if (str(text[i]) == 'Nt'):
                    self.Nt = int(str(text[i+1]))
                if (str(text[i]) == 'Ncolor'):
                    self.Ncolor = int(str(text[i+1]))
                if (str(text[i]) == 'Tpulse'):
                    self.Tpulse = float(str(text[i+1].split()[0]))
                if (str(text[i]) == 'omegac'):
                    self.omegac = float(str(text[i+1].split()[0]))
                if (str(text[i]) == 'phi_CEP'):
                    self.phi_CEP = float(str(text[i+1].split()[0]))
                if (str(text[i]) == 'nenvelope'):
                    self.nenvelope = int(str(text[i+1].split()[0]))
                if (str(text[i]) == 'E0'):
                    self.E0 = float(str(text[i+1].split()[0]))
            if (self.Ncolor >= 2):
                self.omegac = self.omegac*np.ones([self.Ncolor],dtype='float64')
                self.phi_CEP = self.phi_CEP*np.ones([self.Ncolor],dtype='float64')
                self.Tpulse = self.Tpulse*np.ones([self.Ncolor],dtype='float64')
                self.nenvelope = self.nenvelope*np.ones([self.Ncolor],dtype='int32')
                self.E0 = self.E0*np.ones([self.Ncolor],dtype='float64')
                for i in range(Nlen):
                    if (str(text[i]) == 'Tpulse'):
                        temp = text[i+1].split()
                        self.Tpulse = np.array(temp,dtype='float64')
                        if (len(self.Tpulse) != self.Ncolor):
                            print('Error: Number of argmeunt in Tpulse is wrong.')
                            sys.exit()
                    if (str(text[i]) == 'omegac'):
                        temp = text[i+1].split()
                        self.omegac = np.array(temp,dtype='float64')
                        if (len(self.omegac) != self.Ncolor):
                            print('Error: Number of argmeunt in omegac is wrong.')
                            sys.exit()
                    if (str(text[i]) == 'phi_CEP'):
                        temp = text[i+1].split()
                        self.phi_CEP = np.array(temp,dtype='float64')
                        if (len(self.phi_CEP) != self.Ncolor):
                            print('Error: Number of argmeunt in phi_CEP is wrong.')
                            sys.exit()
                    if (str(text[i]) == 'nenvelope'):
                        temp = text[i+1].split()
                        self.nenvelope = np.array(temp,dtype='int32')
                        if (len(self.nenvelope) != self.Ncolor):
                            print('Error: Number of argmeunt in nenvelope is wrong.')
                            sys.exit()
                    if (str(text[i]) == 'E0'):
                        temp = text[i+1].split()
                        self.E0 = np.array(temp,dtype='float64')
                        print(self.E0)
                        if (len(self.E0) != self.Ncolor):
                            print('Error: Number of argmeunt in E0 is wrong.')
                            sys.exit()
        else:
            print('Error: Number of argmeunt is wrong.')
            sys.exit()
        
        print('#=====Print the parmeters')
        print('# cluster_mode =', self.cluster_mode)
        print('# PC_option =', self.PC_option)
        print('# minimal_output =', self.minimal_output)
        print('# a =', self.a, ' [a.u.] =',self.a*aB, ' [nm]')
        print('# Delta =', self.Delta, ' [a.u.] =', self.Delta*Hartree, ' [eV]') 
        print('# dt =', self.dt, ' [a.u.] =', self.dt*Atomtime, ' [fs]') 
        print('# Nt =', self.Nt) 
        print('# Nt*dt =', self.Nt*self.dt, '[a.u.] =', self.Nt*self.dt*Atomtime, '[fs]')
        print('# Number of color Ncolor = ', self.Ncolor)
        if (self.Ncolor == 1):
            print('# Tpulse =', self.Tpulse, ' [a.u.] =', self.Tpulse*Atomtime, ' [fs]') 
            print('# nenvelope =', self.nenvelope) 
            print('# omegac =', self.omegac, ' [a.u.] =', self.omegac*Hartree, ' [eV]') 
            print('# phi_CEP =', self.phi_CEP, ' [2 pi]')
            self.phi_CEP = tpi*self.phi_CEP
            print('# E0 =', self.E0, ' [a.u.] =', self.E0*Atomfield, ' [V/nm]') 
            print('# e*a*E0 =', self.a*self.E0, ' [a.u.] =', self.a*self.E0*Hartree, ' [eV]')
        else:
            for icolor in range(self.Ncolor):
                print('# =====', icolor,'th color ========')
                print('# Tpulse =', self.Tpulse[icolor], ' [a.u.] =', self.Tpulse[icolor]*Atomtime, ' [fs]') 
                print('# nenvelope =', self.nenvelope[icolor]) 
                print('# omegac =', self.omegac[icolor], ' [a.u.] =', self.omegac[icolor]*Hartree, ' [eV]') 
                print('# phi_CEP =', self.phi_CEP[icolor], ' [2 pi]')
                self.phi_CEP[icolor] = tpi*self.phi_CEP[icolor]
                print('# E0 =', self.E0[icolor], ' [a.u.] =', self.E0[icolor]*Atomfield, ' [V/nm]') 
                print('# e*a*E0 =', self.a*self.E0[icolor], ' [a.u.] =', self.a*self.E0[icolor]*Hartree, ' [eV]')

