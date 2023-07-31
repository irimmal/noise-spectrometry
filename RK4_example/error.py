import numpy as np
import matplotlib
matplotlib.use('Agg') # the 'TkAgg' backend broke somehow. not solved yet...
from matplotlib import pyplot as plt
from Ex2_RK4_variables import t_A, E_A, E1, E2, V, dt, T, psi_0
from RK4 import rk4_func
import pandas as pd
import os
import sys
import h5py
 
'''
error analysis for RK4
'''
try:
    save_ERR = bool(sys.argv[1])
except:
    save_ERR = False

H = np.array([[E1, V], [V, E2]], dtype=np.complex64)

N = int(T/dt)
psi = np.zeros(shape=(2, N), dtype=np.complex64, order='F') # argument "order='F'" is important!!!
psi[:, 0] = psi_0

'''calculation: 4th Runge-Kutta'''
### f2py
rk4_func.calculate_rk4(psi, psi_0, H, dt, N)

'''exact'''
t_axis = np.linspace(0, (N-1)*dt, N)

# basis
if not (V==0 and 0.5*((E1-E2) + ((E1-E2)**2 + 4*V**2)**0.5) == 0):
    psi_1 = np.array([-V, 0.5*((E1-E2) + ((E1-E2)**2 + 4*V**2)**0.5)], dtype=np.complex64)
    psi_2 = np.array([0.5*((E1-E2) + ((E1-E2)**2 + 4*V**2)**0.5), V], dtype=np.complex64)
else:
    psi_1 = np.array([0.5*((E1-E2) - ((E1-E2)**2 + 4*V**2)**0.5), V], dtype=np.complex64)
    psi_2 = np.array([-V, 0.5*((E1-E2) - ((E1-E2)**2 + 4*V**2)**0.5)], dtype=np.complex64)

# normalization
psi_1 = psi_1 / (np.conj(psi_1) @ psi_1)**0.5
psi_2 = psi_2 / (np.conj(psi_2) @ psi_2)**0.5

# projection of initial value onto eigenvectors
C_01 = np.conj(psi_0) @ psi_1
C_02 = np.conj(psi_0) @ psi_2

E_eig_1 = 0.5*((E1+E2) - ((E1-E2)**2 + 4*V**2)**0.5)
E_eig_2 = 0.5*((E1+E2) + ((E1-E2)**2 + 4*V**2)**0.5)
psi_exact_1 = C_01*psi_1[0]*np.exp(1/1j * E_eig_1 * t_axis) + C_02*psi_2[0]*np.exp(1/1j * E_eig_2 * t_axis)
psi_exact_2 = C_01*psi_1[1]*np.exp(1/1j * E_eig_1 * t_axis) + C_02*psi_2[1]*np.exp(1/1j * E_eig_2 * t_axis)

'''error analysis'''
ERR = (np.abs(psi[0,:] - psi_exact_1)**2 + np.abs(psi[1,:] - psi_exact_2)**2)**0.5 / (np.abs(psi_exact_1)**2 + np.abs(psi_exact_2)**2)**0.5
err_fit = np.polyfit(t_axis, ERR, deg=1)
print(dt, err_fit[0], err_fit[1])

if save_ERR == True:
    '''write into HDF5'''
    fileName = 'ERR.h5'
    if os.path.isfile(fileName):
        f = h5py.File(fileName, 'a')
    else:
        f = h5py.File(fileName, 'w')

    num_dataset = len(list(f.keys()))
    dset_name = '{0:04}'.format(num_dataset)
    f.create_dataset(dset_name, data=ERR, dtype = np.float64)
    f[dset_name].attrs['dt'] = dt
    f[dset_name].attrs['fit'] = err_fit
    f.close()

'''visualization'''
title_prop = {'fontsize': 16}
text_prop = {'fontsize': 14}

# conversion of units
t_axis_SI = t_A*t_axis

fig = plt.figure()

plt.plot(t_axis_SI, ERR, label=r'error $\Delta C_2(t)$', lw=1, color='r')
plt.plot(t_axis_SI, np.poly1d(err_fit)(t_axis), 'k--')

title_prop = {'fontsize': 16}
text_prop = {'fontsize': 14}
plt.xlabel('time (fs)', text_prop)
plt.ylabel('relative deviation from exact sol.', text_prop)
plt.title(r'($E_1$, $E_2$, V)' + f'={round(E1*1e3/E_A, 2), round(E2*1e3/E_A, 2), round(V*1e3/E_A, 2)} meV', title_prop)
plt.legend()

plt.savefig(r'ERR_E1_E2_V' + f'={round(E1*1e3/E_A, 2), round(E2*1e3/E_A, 2), round(V*1e3/E_A, 2), dt}_meV_cos2.png')
