import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from Ex2_RK4_variables import t_A, E_A, E1, E2, V, dt, T, psi_0
# from RK4 import RK4_func
from RK4 import rk4_func

'''
Schrödinger equation for 2 level system
'''

H = np.array([[E1, V], [V, E2]], dtype=np.complex64)

N = int(T/dt)
psi = np.zeros(shape=(2, N), dtype=np.complex64, order='F') # argument "order='F'" is important!!!
psi[:, 0] = psi_0

'''calculation: 4th Runge-Kutta'''
# ### f2py
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

# probably useful when the Hilbert space is larger
C_0 = np.tensordot(np.conj(psi_0), np.array([psi_1, psi_2]), axes=([0],[1]))
psi_exact_test_tensor = np.tensordot(
                    np.array([np.kron(psi_1, np.exp(1/1j * E_eig_1 * t_axis)).reshape(2, N), 
                            np.kron(psi_2, np.exp(1/1j * E_eig_2 * t_axis)).reshape(2, N)]), 
                    C_0, axes=([0],[0]))

'''visualization'''
title_prop = {'fontsize': 16}
text_prop = {'fontsize': 14}

# conversion of units
t_axis_SI = t_A*t_axis

fig = plt.figure()

# RK4 calculation
plt.plot(t_axis_SI, np.abs(psi[0, :])**2, 'bo', ms=2, label=r'RK4 $|C_1(t)|^2$')
plt.plot(t_axis_SI, np.abs(psi[1, :])**2, 'g^', ms=2, label=r'RK4 $|C_2(t)|^2$')

# exact
plt.plot(t_axis_SI, np.abs(psi_exact_1)**2, label=r'exact $|C_1(t)|^2$', lw=1, color='c')
plt.plot(t_axis_SI, np.abs(psi_exact_2)**2, label=r'exact $|C_2(t)|^2$', lw=1, color='r')

# plt.plot(t_axis, np.abs(psi_exact_test_tensor[0,:])**2, label='exact 0')
# plt.plot(t_axis, np.abs(psi_exact_test_tensor[1,:])**2, label='tensor 1')

# # theoretical approximation, E1 = E2
# import matplotlib.colors as mcolors
# plt.plot(t_axis_SI, np.cos(V*t_axis)**2, label=r'$cos^2(V/\hbar t)$', lw=1, color=mcolors.CSS4_COLORS['lawngreen'])
# plt.plot(t_axis_SI, np.sin(V*t_axis)**2, label=r'$sin^2(V/\hbar t)$', lw=1, color=mcolors.CSS4_COLORS['darkviolet'])


# # theoretical approximation, |E1 - E2| >> V
# plt.plot(t_axis_SI, 1 - 2*V**2/(E1-E2)**2*(1 - np.cos((E1-E2 + 2*V**2/(E1-E2))*t_axis)), label=r'theory (approx)', lw=1, color=mcolors.CSS4_COLORS['lawngreen'])
# plt.plot(t_axis_SI, 2*V**2/(E1-E2)**2*(1 - np.cos((E1-E2 + 2*V**2/(E1-E2))*t_axis)), label=r'theory (approx)', lw=1, color=mcolors.CSS4_COLORS['darkviolet'])


title_prop = {'fontsize': 16}
text_prop = {'fontsize': 14}
plt.xlabel('time (fs)', text_prop)
plt.ylabel('occupation of levels ($|C_1(t)|^2, |C_2(t)|^2$)', text_prop)
plt.title(r'($E_1$, $E_2$, V)' + f'={round(E1*1e3/E_A, 2), round(E2*1e3/E_A, 2), round(V*1e3/E_A, 2)} meV', title_prop)
plt.legend()

plt.savefig(r'occupation_E1_E2_V' + f'={round(E1*1e3/E_A, 2), round(E2*1e3/E_A, 2), round(V*1e3/E_A, 2)}_meV_cos2.png')
# plt.show()
