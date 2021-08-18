This repository contains the software I developed for my thesis for my MSc in Astronomy and Astrophysics.

It requires [MESA](http://mesa.sourceforge.net/) and [GYRE](https://gyre.readthedocs.io/) for creating models and uses [`mesa_reader`](https://github.com/wmwolf/py_mesa_reader) and `h5py` for reading these models.
`numpy, scipy` and `astropy` are used in most of the scripts.
`statsmodels` is used for the inversion routines.
The minimum python version required to run these scripts is Python 3.8.
The code is licensed under the GNU Affero General Public License version 3.

It contains the following tools:

* Running MESA and GYRE simulations: `mesa.py` and `gyre.py` with template inlist files in `inlist_templates`. Requires a choice of template, the parameters for that template (mass, overshooting length scale, ... for MESA; l and m values, model and profile for GYRE) and an optional name for the model. The output of the simulations will be written to the directory `models`. GYRE output can be found in the `gyre` subfolder of the MESA model.  Before creating a MESA model, you should that a MESA work directory has been initialized under `models/workdir/MESA`.
* Computing oscillation kernels with variables (ρ,c) and (N²,c): `kernels.py`. The kernels are written to the HDF files for each of the modes.
* Computing inverted profiles with regularized least squares: `inversion.py` for computing an inversion profile, comparing frequencies from two different models. `inversion_cross_validation.py` for using cross-validation to find the best regularization parameter.
* Application of the inversion routines to the start KIC 10526294 based on models by [Moravveji et al. (2015)](https://doi.org/10.1051/0004-6361/201425290) can be found in `moravveji.py`.
* Some of the functionality has been extracted in other scripts: `datastructures.py` and `mesa_reader_wrapper.py` contain code to interact with MESA and GYRE output. `inversionlib.py` contains the main inversion routines.
* The extra meshing code for the overshoot zone, adapted from [Michielsen et al. (2021)](https://doi.org/10.1051/0004-6361/202039926), as well as the routines for saving models at exact Xc can be found in `run_star_extras.f90`.
