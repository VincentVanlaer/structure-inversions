from math import pi
from inversionlib import invert, cross_validate
from scipy.interpolate import interp1d
from scipy.integrate import cumulative_trapezoid
from datastructures import GyreDescriptor, centered_average
from mesa_reader_wrapper import constants
import numpy as np
import matplotlib.pyplot as plt

g = GyreDescriptor('KIC10526294/gyre_moravveji')

N2_kernels = True
do_cross_validation = False
mu = 1e-5 # overwritten when do_cross_validation is True
x_steps = np.linspace(0, 1, 200)

freq = np.array([
    0.472220,
    0.486192,
    0.500926,
    0.517303,
    0.533426,
    0.552608,
    0.571964,
    0.593598,
    0.615472,
    0.641202,
    0.670600,
    0.701246,
    0.734708,
    0.772399,
    0.812940,
    0.856351,
    0.902834,
    0.954107,
    1.013415
    ])

G = constants.G

end_freqs = {k: v* 2 * pi / 86400 for (k,v) in zip(range(14,33), reversed(freq))}

start_modes = {}

for mode in g.modes:
    mode = mode.load_raw()

    start_modes[mode.attrs['n_g']] = mode

mode = start_modes[1]

R = mode.attrs['R_star']
M = mode.attrs['M_star']
x = mode['x'][:]
rho = mode['rho'][:]
M_r = mode['M_r'][:]
with np.errstate(divide='ignore', invalid='ignore'):
    g = G * M_r * (x * R)**-2
g[0] = 0
P = mode['P'][:]
Gamma_1 = mode['Gamma_1'][:]
rho = mode['rho']
c_sqrd = P / rho * Gamma_1
brunt_N2 = mode['As'][:] / mode['c_1'] * G * M / R**3
freqs = {}

for n_g, mode in start_modes.items():
    freqs[n_g] = mode.attrs['freq'][0] * 2 * pi / 86400

shared_modes = list(start_modes.keys() & end_freqs.keys())[:40]

start_modes = [start_modes[mode] for mode in shared_modes]
start_freqs = np.fromiter((freqs[mode] for mode in shared_modes), dtype=float)
end_freqs = np.fromiter((end_freqs[mode] for mode in shared_modes), dtype=float)

N = len(x_steps) - 1

if N2_kernels:
    kernels = ['kernels/N2_c/N2']
    xname = ['q'] + [f"N{x}" for x in range(N)]
else:
    kernels = ['kernels/rho_c/c', 'kernels/rho_c/rho']
    xname = ['q'] + [f"c{x}" for x in range(N)] + [f"rho{x}" for x in range(N)]

if do_cross_validation:
    errors = np.empty(10)
    for i, mu in enumerate(10**np.linspace(0, -9, 10)):
        errors[i] = cross_validate(start_modes, end_freqs, kernels, x_steps, mu)
        print(errors[i], mu)

    lowest = np.argmin(errors)
    errors = np.empty(20)
    exponent = np.linspace(-lowest + 1, -lowest - 1, 20)

    for i, mu in enumerate(10**exponent):
        errors[i] = cross_validate(start_modes, end_freqs, kernels, x_steps, mu)
        print(errors[i], mu)

    lowest = np.argmin(errors)
    mu = 10**exponent[lowest]

rlsFit = invert(start_modes, end_freqs, kernels, x_steps, mu)

print(rlsFit.summary(xname=xname))

if N2_kernels:
    delta_N2 = np.append(rlsFit.params[1:N+1], rlsFit.params[N])
else:
    delta_c_estimated = np.append(rlsFit.params[1:N+1], rlsFit.params[N])
    delta_rho_estimated = np.append(rlsFit.params[N+1:], rlsFit.params[-1])

predicted_freqs = np.fromiter((freqs[n_g] * (1 + rlsFit.fittedvalues[shared_modes.index(n_g)]) for n_g in shared_modes), dtype=float)

x_axis = 1/start_freqs * 2 * pi / 86400
x_label = 'Angular period ($d^{-1}$)'

plt.plot(x_axis, (end_freqs-start_freqs)/start_freqs, '-o', markersize=3, label='Target')
plt.plot(x_axis, (predicted_freqs-start_freqs)/start_freqs, '-o', markersize=3, label='Predicted')
plt.plot(x_axis, (end_freqs-predicted_freqs)/start_freqs, '-o', markersize=3, label='Predicted')
plt.axhline(0, linestyle='--', label='Starting model')
plt.xlabel(x_label)
plt.ylabel('Relative frequency difference')
plt.legend()
plt.show()

if N2_kernels:
    delta_N2 = interp1d(x_steps[:-1], delta_N2[:-1], kind='previous', fill_value='extrapolate')(x)
    plt.plot(x, delta_N2)
    plt.plot(x, brunt_N2 * R**3 / G / M)
    plt.plot(x, brunt_N2 * R**3 / G / M + delta_N2)
    plt.show()

else:
    plt.plot(x_steps, delta_rho_estimated, drawstyle='steps-post', label='Inverted density')
    plt.plot(x_steps, delta_c_estimated, drawstyle='steps-post', label='Inverted sound speed')
    plt.legend()
    plt.show()


    delta_der_rho = centered_average(x_steps, delta_rho_estimated)
    delta_c_estimated = interp1d(x_steps[:-1], delta_c_estimated[:-1], kind='previous', fill_value='extrapolate')
    delta_rho_estimated = interp1d(x_steps[:-1], delta_rho_estimated[:-1], kind='previous', fill_value='extrapolate')


    delta_der_rho = interp1d(x_steps[:-1], delta_der_rho[:-1], kind='previous', fill_value='extrapolate')(x)
    g = g * R**2 / G / M
    c_sqrd = c_sqrd* G * M / R
    delta_g = 4 * pi * G * x**-2 * cumulative_trapezoid(x**2 * delta_rho_estimated(x) * rho * R**3 / M, x, initial=0)
    delta_N2 = delta_g / g * (brunt_N2 - g**2 / c_sqrd) - 2 * g**2 / c_sqrd * delta_c_estimated(x) - g * delta_der_rho

    plt.plot(x, delta_N2)
    plt.plot(x, brunt_N2 * R**3 / G / M)
    plt.show()
