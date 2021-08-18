from __future__ import annotations
from functools import cached_property
from typing import List, Dict
from mesa_reader_wrapper import MesaLog
from pathlib import Path
import astropy.table
import h5py
import numpy as np

model_path = Path('models')


def centered_average(x,y):
    h = x[1:] - x[:-1]
    s = (y[1:] - y[:-1]) / h

    dx = np.empty(x.shape)
    dx[0] = s[0]
    dx[1:-1] = (s[:-1]*h[1:] + s[1:]*h[:-1]) / (h[:-1] + h[1:])
    dx[-1] = s[-1]

    return dx


def set_global_model_path(new_model_path: Path) -> None:
    global model_path
    model_path = new_model_path


class MesaModelDescriptor:
    def __init__(self, identifier: str) -> None:
        self.identifier = identifier

    @cached_property
    def profiles(self) -> Dict[MesaProfileDescriptor]:
        profiles = {}

        for f in model_path.joinpath(self.identifier, 'LOGS').glob('profile*.data'):
            index = int(f.stem[len('profile'):])
            profiles[index] = MesaProfileDescriptor(self, index)

        return profiles


class MesaProfileDescriptor:
    def __init__(self, model: MesaModelDescriptor, index: int) -> None:
        self.model = model
        self.index = index

    @property
    def identifier(self) -> str:
        return f'{self.model.identifier}/LOGS/profile{self.index}.data'

    def load(self) -> MesaLog:
        return MesaLog(model_path.joinpath(self.identifier))


def _load_gyre_as_table(f: h5py.File) -> astropy.table.Table:
    cols = {}

    for key, val in f.items():
        if not isinstance(val, h5py.Group):
            cols[key] = val[...]
    tab = astropy.table.Table(cols, meta=dict(zip(f.attrs.keys(), f.attrs.values())))

    complex_dtype = np.dtype([('re', '<f8'), ('im', '<f8')])

    for colname in tab.colnames:

        if tab[colname].dtype == complex_dtype:
            tab[colname] = tab[colname]['re'] + 1j*tab[colname]['im']

    for atrname in tab.meta:

        if tab.meta[atrname].dtype == complex_dtype:
            tab.meta[atrname] = tab.meta[atrname]['re'] + 1j*tab.meta[atrname]['im']

    return tab


class GyreDescriptor:

    def __init__(self, identifier: str) -> None:
        self.identifier = identifier

    @classmethod
    def from_mesa_profile(cls, profile: MesaProfileDescriptor, name: str):
        return cls(f'{profile.model.identifier}/gyre/{profile.index}/{name}')

    def load_summary(self) -> astropy.table.Table:
        return _load_gyre_as_table(h5py.File(model_path.joinpath(self.identifier, 'summary_ad_g_modes.HDF')))

    @cached_property
    def modes(self) -> List[GyreOscillationMode]:
        return [GyreOscillationMode(self, f.name)
                for f in sorted(model_path.joinpath(self.identifier)
                .glob('ad_*.HDF'))]


class GyreOscillationMode:

    def __init__(self, gyre: GyreDescriptor, identifier: str):
        self.gyre = gyre
        self.identifier = identifier

    def load(self) -> astropy.table.Table:
        with self.load_raw() as f:
            return _load_gyre_as_table(f)

    def load_raw(self) -> h5py.File:
        return h5py.File(model_path.joinpath(self.gyre.identifier, self.identifier), 'r+')
