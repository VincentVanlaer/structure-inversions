from __future__ import annotations
from typing import List
import statsmodels.api as sm
import numpy as np
import h5py
from math import pi
from scipy.interpolate import interp1d
from scipy.integrate import trapezoid, cumulative_trapezoid
from scipy.linalg import block_diag
from mesa_reader_wrapper import constants

G = constants.G

def prepare_kernels(modes: List[h5py.File], kernels: List[str]):

    prepared_kernels = []
    mode_grids = []
    freqs = np.empty(len(modes))

    for i, mode in enumerate(modes):
        freqs[i] = mode.attrs['freq'][0] * 2 * pi / 86400
        mode_grids.append(mode['x'][:])
        prepared_kernels.append((*(mode[kernel] for kernel in kernels), ))

    return prepared_kernels, mode_grids, freqs


def rls(observations, mode_grids, constant, kernels, x_steps, mu):
    n_modes = len(kernels)
    n_variables = len(kernels[0])
    n_segments = len(x_steps) - 1

    G = np.empty((n_modes, n_segments * n_variables))

    for j, kernel_set in enumerate(kernels):
        x = mode_grids[j]
        xn = np.searchsorted(x, x_steps, side='right') - 1
        xn[0] = 0
        for k in range(n_variables):
            G[j, k * n_segments:(k + 1) * n_segments] = np.diff(kernel_set[k][xn])

    L = np.diag(np.full(n_segments,-2)) + np.diag(np.ones(n_segments-1),1) + np.diag(np.ones(n_segments-1),-1)
    L[0,0] = L[0,2] = L[-1,-1] = L[-1,-3] = 1
    L[0,1] = L[-1,-2] = -2

    L = L * n_segments * n_modes  # balance number of smoothing parameters and observations

    L = block_diag(*(L for _ in range(n_variables)))

    X = np.concatenate((G, mu * L))
    X = np.hstack((np.concatenate((constant, np.zeros(n_segments * n_variables)))[:,np.newaxis], X))

    return sm.OLS(np.concatenate((observations, np.zeros(n_segments * n_variables))), X, hasconst=False).fit()


def averaging_kernels(moore_penrose, kernels):
    kernels = np.concatenate((np.array([kernel[1] for kernel in kernels]), np.array([kernel[2] for kernel in kernels])))
    N = len(kernels)

    return moore_penrose[:,:N] @ np.array(kernels)


def add_conservation_law(modes: List[h5py.File], kernels: List[str]):
    conserved = []
    mode = modes[0]
    grid = mode['x'][:]
    R_star = mode.attrs['R_star']
    M_star = mode.attrs['M_star']
    rho = mode['rho']

    for kernel in kernels:
        if kernel.endswith('rho'):
            conserved.append(cumulative_trapezoid(grid**2*rho / M_star * R_star**3, grid, initial=0))
        else:
            conserved.append(np.zeros(grid.shape))

    return [(grid, conserved, 0)]


def invert(modes: List[h5py.File], observed_freqs: np.array, kernels: List[str], x_steps, mu):
    kernels_prepared, mode_grids, model_freqs = prepare_kernels(modes, kernels)

    delta = observed_freqs / model_freqs - 1

    grid, kernels_conserved, delta_conserved = zip(*add_conservation_law(modes, kernels))
    n_conserved = len(delta_conserved)
    delta = np.append(delta, delta_conserved)
    kernels_prepared.extend(kernels_conserved)
    mode_grids = mode_grids.copy()
    mode_grids.extend(grid)

    integrated_kernels = [(*(cumulative_trapezoid(kernel, grid, initial=0) for kernel in kernel_set),) for grid, kernel_set in zip(mode_grids, kernels_prepared)]
    constant = np.array([1]*len(modes) + [0]*n_conserved)

    return rls(delta, mode_grids, constant, integrated_kernels, x_steps, mu)


def cross_validate(modes: List[h5py.File], observed_freqs: np.array, kernels: List[str], x_steps, mu):
    kernels_prepared, mode_grids, model_freqs = prepare_kernels(modes, kernels)

    delta = observed_freqs / model_freqs - 1

    grid, kernels_conserved, delta_conserved = zip(*add_conservation_law(modes, kernels))
    n_conserved = len(delta_conserved)
    delta = np.append(delta, delta_conserved)
    kernels_prepared.extend(kernels_conserved)
    mode_grids = mode_grids.copy()
    mode_grids.extend(grid)

    n_variables = len(kernels)
    n_modes = len(modes)
    n_segments = len(x_steps) - 1

    integrated_kernels = [(*(cumulative_trapezoid(kernel, grid, initial=0) for kernel in kernel_set),) for grid, kernel_set in zip(mode_grids, kernels_prepared)]
    constant = np.array([1]*len(modes) + [0]*n_conserved)

    error = 0

    for test_mode in range(n_modes):
        rlsFit = rls(np.delete(delta, test_mode),
                mode_grids[:test_mode] + mode_grids[test_mode+1:],
                np.delete(constant, test_mode),
                integrated_kernels[:test_mode] + integrated_kernels[test_mode + 1:],
                x_steps, mu)
        q = rlsFit.params[0]

        predicted = q

        for variable in range(n_variables):
            inverted = rlsFit.params[1 + n_segments * variable:1 + n_segments * (variable + 1)]
            inverted = interp1d(x_steps[:-1], inverted, kind='previous', fill_value='extrapolate')

            predicted += trapezoid(inverted(mode_grids[test_mode]) * kernels_prepared[test_mode][variable], mode_grids[test_mode])

        error += (predicted - (observed_freqs[test_mode] / model_freqs[test_mode] - 1))**2

    return error


def discretization_only_constant(n):
    yield ('Constant', np.linspace(0, 1, n + 1))

def discretization_constant_focussed(n, mixed_zone_start, mixed_zone_end):
    yield ('Constant', np.linspace(0, 1, n + 1))

    zone = np.linspace(mixed_zone_start*0.99, mixed_zone_end*1.01, int(0.9*n))
    zone = np.concatenate((
        zone,
        np.linspace(0, 1, int(0.1 * n) + 1)))
    zone.sort()

    yield ('Focussed', zone)

    zone = np.linspace(0.2, 0.3, int(0.8*n))
    zone = np.concatenate((
        zone,
        np.linspace(0, 1, int(0.2 * n) + 1)))
    zone.sort()

    yield ('Concentrated #1', zone)

    zone = np.linspace(0.5, 0.6, int(0.8*n))
    zone = np.concatenate((
        zone,
        np.linspace(0, 1, int(0.2 * n) + 1)))
    zone.sort()

    yield ('Concentrated #2', zone)

def discretization_only_focussed(n, mixed_zone_start, mixed_zone_end):
    zone = np.linspace(mixed_zone_start*0.99, mixed_zone_end*1.01, int(0.9*n))
    zone = np.concatenate((
        zone,
        np.linspace(0, 1, int(0.1 * n) + 1)))
    zone.sort()

    yield ('Focussed', zone)
