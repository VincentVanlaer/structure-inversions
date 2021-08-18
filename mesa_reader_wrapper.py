from functools import cached_property
from scipy import interpolate
from mesa_reader import MesaData as UnwrappedMesaData
from collections.abc import Iterable
from pathlib import Path


# as defined in MESA
class constants:
    G = 6.67430e-8  # cm^3 / g / s^2
    M_sun = 1.3271244e26 / G  # g
    R_sun = 6.957e10  # cm


class MesaLog(UnwrappedMesaData):

    def __init__(self, log_path: Path, units: str = 'cgs'):
        super().__init__(log_path, 'log')

        self.unit_converter = MesaLog.__cgs_convert if units == 'cgs' else lambda x: x

    @staticmethod
    def __cgs_convert(name: str, array):
        try:
            return {
                'rmid': lambda x: x * constants.R_sun,
                'radius': lambda x: x * constants.R_sun,
                'mass': lambda x: x * constants.M_sun,
                'mmid': lambda x: x * constants.M_sun,
            }[name](array)
        except KeyError:
            return array

    def _log_version(self, key):
        log_prefixes = ['log_', 'log', 'lg_', 'lg']
        for prefix in log_prefixes:
            if self.in_data(prefix + key):
                return prefix + key
            if self.in_data(prefix + key.capitalize()):
                return prefix + key.capitalize()

    def __getattr__(self, x):
        return self.unit_converter(x, self.data(x))

    @cached_property
    def R(self):
        return self.unit_converter('radius', self.data('radius')[0])

    @cached_property
    def M(self):
        return self.unit_converter('mass', self.header('initial_mass'))

    @cached_property
    def grid(self):
        return self.rmid / self.R


class MesaLogRegridded:

    def __init__(self, original: MesaLog, grid):
        self.original = original
        self.grid = grid

    def __interpolate(self, data):
        return interpolate.interp1d(self.original.grid, data, fill_value='extrapolate')(self.grid)

    def __getattr__(self, x):
        try:
            data = getattr(self.original, x)
            if isinstance(data, Iterable):
                return self.__interpolate(data)
            else:
                return data
        except KeyError:
            return getattr(self.original, x)
