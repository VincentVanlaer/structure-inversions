import json
from argparse import ArgumentParser
from pathlib import Path
from math import pi
import numpy as np
from mesa_reader_wrapper import constants
from datastructures import MesaModelDescriptor, GyreDescriptor, set_global_model_path, model_path
from inversionlib import cross_validate, discretization_constant_focussed

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

start_model_name, start_profile_index = args.start_model.rsplit('/', maxsplit=1)
start_profile = MesaModelDescriptor(start_model_name).profiles[int(start_profile_index)]
start_gyre = GyreDescriptor.from_mesa_profile(start_profile, args.gyre)
start = start_profile.load()

end_model_name, end_profile_index = args.end_model.rsplit('/', maxsplit=1)
end_profile = MesaModelDescriptor(end_model_name).profiles[int(end_profile_index)]
end_gyre = GyreDescriptor.from_mesa_profile(end_profile, args.gyre)

start_modes = {}

for mode in start_gyre.modes:
    mode = mode.load_raw()

    start_modes[mode.attrs['n_g']] = mode

end_modes = end_gyre.load_summary()
end_freqs = {}

for n_g, freq in zip(end_modes['n_g'], end_modes['freq']):
    end_freqs[n_g] = freq.real * 2 * pi / 86400

shared_modes = list(start_modes.keys() & end_freqs.keys())

start_modes = [start_modes[mode] for mode in shared_modes]
end_freqs = [end_freqs[mode] for mode in shared_modes]


mixed_zone_indices = np.logical_and(start.brunt_N2 > 1e-6, (start.rmid / start.R) < 0.5)
mixed_zone_start = mixed_zone_end = None

for index, elem in enumerate(mixed_zone_indices[::-1]):
    if elem and not mixed_zone_start:
        mixed_zone_start = start.rmid[-index - 1] / start.R
    if not elem and mixed_zone_start:
        mixed_zone_end = start.rmid[-index - 1] / start.R
        break

print(f"Identified mixed region as {mixed_zone_start} to {mixed_zone_end}")

all_distributions = []

for iterator in (discretization_constant_focussed(100, mixed_zone_start, mixed_zone_end),):
    distributions = {}
    for name, distribution in iterator:
        errors = np.empty(8)
        for i, mu in enumerate(10**np.linspace(0, -7, 8)):
            errors[i] = cross_validate(start_modes, end_freqs, ['kernels/rho_c/rho', 'kernels/rho_c/c'], distribution, mu)
            print(errors[i], mu)

        lowest = np.argmin(errors)
        errors = np.empty(20)
        exponent = np.linspace(-lowest + 1, -lowest - 1, 20)

        for i, mu in enumerate(10**exponent):
            errors[i] = cross_validate(start_modes, end_freqs, ['kernels/rho_c/rho', 'kernels/rho_c/c'], distribution, mu)
            print(errors[i], mu)

        lowest = np.argmin(errors)
        distributions[name] = (errors[lowest], 10**exponent[lowest])

    all_distributions.append(distributions)

print(all_distributions)
