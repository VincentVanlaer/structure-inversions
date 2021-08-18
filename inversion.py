from argparse import ArgumentParser
from pathlib import Path
from math import pi
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from mesa_reader_wrapper import MesaLogRegridded, constants
from datastructures import MesaModelDescriptor, GyreDescriptor, set_global_model_path, model_path
from inversionlib import invert, discretization_only_focussed

parser = ArgumentParser()

parser.add_argument('--model-dir',
                    type=Path,
                    default=model_path)
# where the gyre output can be found, relative to the mesa model
parser.add_argument('--gyre', type=str, required=True)
# mesa .LOG for starting model
parser.add_argument('--start-model', type=str, required=True)
# mesa .LOG for target model
parser.add_argument('--end-model', type=str, required=True)

args = parser.parse_args()

set_global_model_path(args.model_dir)

G = constants.G

model_name, profile_index = args.start_model.rsplit('/', maxsplit=1)
start_profile = MesaModelDescriptor(model_name).profiles[int(profile_index)]
start_gyre = GyreDescriptor.from_mesa_profile(start_profile, args.gyre)
start = start_profile.load()

model_name, profile_index = args.end_model.rsplit('/', maxsplit=1)
end_profile = MesaModelDescriptor(model_name).profiles[int(profile_index)]
end_gyre = GyreDescriptor.from_mesa_profile(end_profile, args.gyre)
end = MesaLogRegridded(end_profile.load(), start.grid)

brunt_N2_start = start.brunt_N2 * start.R**3 / start.M / G
brunt_N2_end = end.brunt_N2 * end.R**3 / end.M / G

with np.errstate(invalid='ignore'):
    brunt_N_start = np.sqrt(brunt_N2_start)
    brunt_N_end = np.sqrt(brunt_N2_end)

np.nan_to_num(brunt_N_start, copy=False)
np.nan_to_num(brunt_N_end, copy=False)

start_modes = {}

for mode in start_gyre.modes:
    mode = mode.load_raw()

    start_modes[mode.attrs['n_g']] = mode

end_modes = end_gyre.load_summary()
end_freqs = {}
freqs = {}

for n_g, freq in zip(end_modes['n_g'], end_modes['freq']):
    end_freqs[n_g] = freq.real * 2 * pi / 86400

for n_g, mode in start_modes.items():
    freqs[n_g] = mode.attrs['freq'][0] * 2 * pi / 86400

shared_modes = list(start_modes.keys() & end_freqs.keys())[:40]

start_modes = [start_modes[mode] for mode in shared_modes]
start_freqs = np.fromiter((freqs[mode] for mode in shared_modes), dtype=float)
end_freqs = np.fromiter((end_freqs[mode] for mode in shared_modes), dtype=float)

mixed_zone_indices = np.logical_and(start.brunt_N2 > 1e-6, (start.rmid / start.R) < 0.5)
mixed_zone_start = mixed_zone_end = None

for index, elem in enumerate(mixed_zone_indices[::-1]):
    if elem and not mixed_zone_start:
        mixed_zone_start = start.rmid[-index - 1] / start.R
    if not elem and mixed_zone_start:
        mixed_zone_end = start.rmid[-index - 1] / start.R
        break

x_steps = next(discretization_only_focussed(200, mixed_zone_start, mixed_zone_end))[1]
mu = 1e-5

N = len(x_steps) - 1

rlsFit = invert(start_modes, end_freqs, ['kernels/zero/N2_c/N2'], x_steps, 0.000000002423987193956)

print(rlsFit.summary())

delta_N2_estimated = np.append(rlsFit.params[1:N+1], rlsFit.params[N])

delta_N2_estimated = interp1d(x_steps[:-1], delta_N2_estimated[:-1], kind='previous', fill_value='extrapolate')
delta_N2 = delta_N2_estimated(start.grid)

plt.plot(start.grid, (brunt_N2_end - brunt_N2_start), label='$N^2$')
plt.plot(start.grid, delta_N2, label='Inverted $N^2$')
plt.xlim([0.15,0.22])
fig = plt.gcf()
fig.set_size_inches((fig.get_size_inches()[0], fig.get_size_inches()[1] / 1.5))
fig.supylabel('Dimensionless difference')
fig.supxlabel('Radial coordinate (R$_*$)')
plt.legend()
plt.show()

predicted_freqs_inverted = {}
predicted_freqs = {}

for n_pg in shared_modes:
    predicted_freqs[n_pg] = freqs[n_pg] * (1 + rlsFit.fittedvalues[shared_modes.index(n_pg)])

x_axis = 1/start_freqs
x_label = 'Angular period (dimensionless)'

predicted_freqs = np.fromiter((predicted_freqs[mode] for mode in shared_modes), dtype=float)
plt.plot(x_axis, (end_freqs-start_freqs)/start_freqs, '-o', markersize=3, label='Target model')
plt.plot(x_axis, (predicted_freqs-start_freqs)/start_freqs, '-o', markersize=3, label='Predicted')
plt.axhline(0, linestyle='--', label='Starting model')
plt.xlabel(x_label)
plt.ylabel('Relative frequency difference')
plt.legend()
fig = plt.gcf()
fig.set_size_inches((fig.get_size_inches()[0], fig.get_size_inches()[1] / 1.5))
plt.show()
