from math import pi
from datastructures import centered_average, GyreDescriptor, MesaModelDescriptor, set_global_model_path, model_path
from mesa_reader_wrapper import constants
from scipy.interpolate import interp1d
from pathlib import Path
from argparse import ArgumentParser
import scipy
import warnings
import numpy as np
import scipy.integrate as integrate

warnings.filterwarnings('ignore', message='The following arguments have no effect for a chosen solver')

G = constants.G

def rho_c_kernels(mode):

    # general properties
    l = mode.meta['l']
    R_star = mode.meta['R_star']
    M_star = mode.meta['M_star']
    freq = mode.meta['freq'].real * 2 * pi / 86400

    # profiles
    grid = mode['x']
    x = R_star * grid
    xi_h = R_star * mode['xi_h'].real
    xi_r = R_star * mode['xi_r'].real
    P = mode['P']
    rho = mode['rho']
    Gamma_1 = mode['Gamma_1']
    c_sqrd = P / rho * Gamma_1
    M_r = mode['M_r']
    with np.errstate(divide='ignore', invalid='ignore'):
        g = G * M_r * x**-2
    g[0] = 0
    L = l * (l+1)
    eul_phi = mode['eul_phi'].real * G * M_star / R_star
    deul_phi = mode['deul_phi'].real * G * M_star / R_star**2

    # Derivatives
    dxi_r = centered_average(x, xi_r)
    drho = centered_average(x, rho)
    with np.errstate(divide='ignore', invalid='ignore'):
        chi = centered_average(x, x**2*xi_r)/x**2 - (xi_h * l * (l + 1) / x)
    chi[0] = 0  # first derivative is always zero in center (symmetry)

    S = np.trapz((xi_r**2 + L*xi_h**2)*rho*x**2, grid)  # grid is used here to make the kernels dimensionless

    df_x_up = np.diff(x)
    df_x_down = df_x_up[::-1]

    # integrals for Δrho kernel when ΔdP and Δdrho have been transformed to Δrho

    int_1 = integrate.cumulative_trapezoid((xi_r * chi * rho)[::-1], dx=df_x_down, initial=0)[::-1]
    int_2 = integrate.cumulative_trapezoid((xi_r**2 * drho)[::-1], dx=df_x_down, initial=0)[::-1]

    c_kernel = 2 * rho * c_sqrd * chi**2 * x**2
    rho_kernel = 2 * (- 0.5 * freq**2 * (xi_r**2 + l*(l+1)*xi_h**2)*rho*x**2
        + 0.5*c_sqrd * chi**2 *rho * x**2
        - G * M_r * rho * xi_r * chi - 4 * pi * G * rho * x**2 * int_1
        + G * M_r * rho * xi_r * dxi_r + 2 * G * pi*x**2 * rho**2*xi_r**2
        - 2*pi * G * rho * x**2 * int_2
        + (L * x * xi_h * eul_phi + deul_phi * xi_r * x**2) * rho
            )

    alpha = integrate.trapezoid(x**2*rho*rho_kernel) / integrate.trapezoid(x**2*rho*x**2*rho)
    complementary = alpha*x**2*rho

    rho_kernel[0] = 0
    # dimensionless (except for integration)
    rho_kernel /= (S * freq**2 * 2)
    complementary /= (S * freq**2 * 2)
    c_kernel /= (S * freq**2 * 2)

    return c_kernel, rho_kernel, complementary


def N2_c_kernels(rho_kernel, c_kernel, mode, minimizer=None):
    # general properties
    R_star = mode.meta['R_star']
    M_star = mode.meta['M_star']

    # profiles
    grid = mode['x']
    x = R_star * grid
    P = mode['P']
    rho = mode['rho']
    Gamma_1 = mode['Gamma_1']
    c_sqrd = P / rho * Gamma_1
    M_r = mode['M_r']
    with np.errstate(divide='ignore', invalid='ignore'):
        g = G * M_r * x**-2
    g[0] = 0
    N2 = mode['As'] / mode['c_1'] * G * M_star / R_star**3

    with np.errstate(divide='ignore', invalid='ignore'):
        b_1 = -2 / x * g + 4 * pi * G * rho; b_1[0] = 4 / 3 * pi * G * rho[0]
        b_2 = 4*pi*G*x**2*rho
        b_3 = x**-2 * (N2 / g - g / c_sqrd)

        b_4 = g

        c_1 = b_1 / b_2 ; c_1[0] = c_1[1]
        c_2 = 1 / b_2 ; c_2[0] = c_2[1]
        c_3 = b_3
        c_4 = b_4 / b_2 ; c_4[0] = c_4[1]

        intermediate_2 = c_2*rho_kernel ; intermediate_2[0] = 0

        d_1 = centered_average(x, c_1) - c_3 ; # a_1
        d_2 = centered_average(x, intermediate_2) # a_0
        d_3 = c_1 + centered_average(x, c_4) ; # a_3
        d_4 = c_4 ; # a_4

        d_1 *= R_star**3
        d_2 *= G * M_star
        d_3 *= R_star**2
        d_4 *= R_star

        e_1 = d_2 / d_4 ; e_1[0] = 0; e_1[1] = e_1[2] / 2
        e_2 = -d_1 / d_4 ; e_2[0] = e_2[1]
        e_3 = -d_3 / d_4 ; e_3[0] = e_3[1]

    if not np.isfinite(e_1).all():
        raise FloatingPointError('Invalid value in e_1')
    if not np.isfinite(e_2).all():
        raise FloatingPointError('Invalid value in e_2')
    if not np.isfinite(e_3).all():
        raise FloatingPointError('Invalid value in e_3')

    e_1 = interp1d(x / R_star, e_1, fill_value='extrapolate')
    e_2 = interp1d(x / R_star, e_2, fill_value='extrapolate')
    e_3 = interp1d(x / R_star, e_3, fill_value='extrapolate')

    def ivp(t, y):
        res = (y[1], e_1(t) + e_2(t) * y[0] + e_3(t) * y[1])
        return res

    def jac(t, y):  # RK45 doesn't use this but other integration methods do
        return ((0, 1), (e_2(t), e_3(t)))

    y = None

    def C_minimizer(y, y_prime):
        y_f = y * R_star**3 / G / M_star
        y_prime_f = y_prime * R_star**2 / G / M_star
        f = y_f * b_3
        f[0] = f[1]
        integral = integrate.trapezoid(f, x)

        return integrate.trapezoid((y_f * b_1 + b_2 * (integral - integrate.cumulative_trapezoid(f, x, initial=0)) + g * y_prime_f + 0*x**2 * rho / M_star * R_star - rho_kernel)**2, x / R_star)

    def to_minimize(y_0):
        nonlocal y
        res = integrate.solve_ivp(ivp, (0,1), (y_0[0], 0), jac=jac, method='RK23', max_step=0.0001)
        y = interp1d(res.t, res.y[0], fill_value='extrapolate')(x / R_star)
        y_prime = interp1d(res.t, res.y[1], fill_value='extrapolate')(x / R_star)

        if minimizer:
            res = minimizer(y, y_prime)
        else:
            res = C_minimizer(y, y_prime)

        return res

    scipy.optimize.minimize(to_minimize, [0], method='Nelder-Mead')

    y2 = c_kernel - 2 * g**2 / c_sqrd * y * R_star**3 / G / M_star

    return y, y2


def write_kernels(gyre: GyreDescriptor, kernel_type: str):

    for mode in gyre.modes:
        mode_astropy = mode.load()
        mode_h5py = mode.load_raw()

        kernels = mode_h5py.require_group('kernels')

        if kernel_type in ['rho_c', 'N2']:
            c_kernel, rho_kernel, complementary = rho_c_kernels(mode_astropy)

            kernels_rho_c = kernels.require_group('rho_c')
            if 'c' in kernels_rho_c:
                del kernels_rho_c['c']
            if 'rho' in kernels_rho_c:
                del kernels_rho_c['rho']
            if 'complementary' in kernels_rho_c:
                del kernels_rho_c['complementary']

            kernels_rho_c.create_dataset('c', data=c_kernel)
            kernels_rho_c.create_dataset('rho', data=rho_kernel)
            kernels_rho_c.create_dataset('complementary', data=rho_kernel)

            if kernel_type == 'N2':
                N2_kernel, c_kernel = N2_c_kernels(rho_kernel, c_kernel, mode_astropy)

                kernels_N2_c = kernels.require_group('N2_c')
                if 'c' in kernels_N2_c:
                    del kernels_N2_c['c']
                if 'N2' in kernels_N2_c:
                    del kernels_N2_c['N2']

                kernels_N2_c.create_dataset('c', data=c_kernel)
                kernels_N2_c.create_dataset('N2', data=N2_kernel)

        mode_h5py.close()


if __name__ == '__main__':
    parser = ArgumentParser()

    parser.add_argument('--model-dir',
                        type=Path,
                        default=model_path)
    # where the gyre output can be found, relative to the mesa model
    parser.add_argument('--gyre', type=str, required=True)
    # mesa .LOG for starting model
    parser.add_argument('--start-model', type=str, required=True)
    # Additional plots
    parser.add_argument('--plot', action='store_true')
    # kernel choice
    parser.add_argument('--kernel-type', choices=['rho_c', 'N2'], default='rho_c')

    args = parser.parse_args()

    set_global_model_path(args.model_dir)

    model_name, profile_index = args.start_model.rsplit('/', maxsplit=1)
    start_profile = MesaModelDescriptor(model_name).profiles[int(profile_index)]
    start_gyre = GyreDescriptor.from_mesa_profile(start_profile, args.gyre)

    write_kernels(start_gyre, args.kernel_type)
