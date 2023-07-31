import numpy as np
'''
Hartree atomic units
e = m_e = hbar = k_e = 1
energy (1 E_A ~ 27.211 eV)
time (t_A ~ 2.419e-17)
'''

'''const'''
t_A = 2.419e-2 # atomic unit -> fs
E_A = 27.211386245988 # atomic unit -> eV


E1 = 1
E2 = 0.5
V = 1 # off-diagonal

# T = 1000 # total simulation time, for error analysis
T = 10 # total simulation time, fo =r plotting
dt = 0.1 # time step

psi_0 = np.array([1, 0], dtype=np.complex64, order='F')
# psi_0 = np.array([0.5, 3**0.5/2], np.complex64)