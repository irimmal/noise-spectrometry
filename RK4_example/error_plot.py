import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
from Ex2_RK4_variables import T, t_A
import h5py

col_name = ['dt', 'slope', 'intercept']
data = [
    pd.read_csv('global_ERR', delimiter=' ', names=col_name)
]

title_prop = {'fontsize': 16}
text_prop = {'fontsize': 14}

'''4th power fit'''
fig = plt.figure(dpi=150)
ROI = (data[0]['dt'] < 0.5)

for i, datum in enumerate(data):
    plt.plot(datum['dt'][ROI], datum['slope'][ROI], 'bo', label='observation')

dt_list = np.linspace(0,0.5,100)

func = lambda x, a: a*x**4
curve_fit_res = curve_fit(func, data[0]['dt'][ROI], data[0]['slope'][ROI], p0=[1e-4])
print(curve_fit_res)
plt.plot(dt_list, func(dt_list, *curve_fit_res[0]), 'r--', label=r'$(\Delta t)^4$ fitting')

plt.xlabel(r'$\Delta t$ (atomic unit)', text_prop)
plt.ylabel('slope (fitting result)', text_prop)
# plt.title('relative global error / T', title_prop)
plt.title('relative global error / T, close-up', title_prop)

plt.legend()
plt.tight_layout()
# plt.savefig('rel_global_error.png')
plt.savefig('rel_global_error_close-up.png')

'''plot typical data and linear fitting of slopes'''

fileName = 'ERR.h5'
f = h5py.File(fileName, 'r')
data = []
data_dt = []
data_fit = []
for name in f.keys():
    data.append(f[name][:])
    data_dt.append(f[name].attrs['dt'])
    data_fit.append(f[name].attrs['fit'])
f.close()


title_prop = {'fontsize': 16}
text_prop = {'fontsize': 14}

fig_ERR = plt.figure(dpi=150)

for i, datum in enumerate(data[:6]):
    t_axis = np.linspace(0, T-data_dt[i], datum.size)
    plt.plot(t_axis, np.poly1d(data_fit[i])(t_axis), 'k--', )
    plt.plot(t_axis, datum, label=r'$\Delta$t=' + f'{data_dt[i]}')

plt.xlabel(f'time (a.u.) = ({round(t_A, 3)} fs)', text_prop)
plt.ylabel('relative deviation from exact sol.', text_prop)
plt.title('error accumulation rate', title_prop)
plt.legend()
plt.savefig('RK_slope_fit_typical.png')
