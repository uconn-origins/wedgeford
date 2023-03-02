# wedgeford
Code repository for the Warm Emission Due to Gas Envelopes Falling Onto Rotating Disks (WEDGEFORD) project

WEDGEFORD is a framework for generating thermochemical models of embedded protostellar objects. The code uses RADMC3D to process three dimensional models of disks and their envelopes in order to calculate the dust temperatures, radiation fields, and gas temperatures and inputs into astrochemical computations. The current [chemistry module](./chemistry) writes inputs for a chemical code based on Fogel et al. (2011).

Libraries used:
radmc3d (fortran executable)
radmcpy (python library)
optool (fortran executable)


## Parameters
The code uses a set of parameters to self-consistently set the information to be used for the physical, thermal, and chemical computations. Sample parameter libraries can be found in [`models/templates/param_library.py`](./models/templates/param_library.py)

To initialize a new object, you just need to pass a dictionary of parameters and an output directory path to the `initialize_model` function. 

To load a model from an existing model directory, the path of the directory can be passed into `load_model`

### Physical Model
The basic physical model is a disk and envelope component with two populations of dust grains (small and large), with only large grains settled in the disk. 

There is the option to load in the density distribution of dust and gas from FARGO data ['models/import_fargo.py'](./models/import_fargo.py)

Further details for how the physical model is generated is in [`models/models.py` ](./models/models.py)

### Radiative Transfer
Radiative transfer of thermal and high energy photons (UV + X-ray) is done with RADMC-3D. How the model interfaces with the RADMC-3D outputs can be found in ['models/outputs.py'](./models/outputs.py).
Models for setting the radiation sources is in ['models/radation.py'](./models/radiation.py) and modules for the dust properties (which can be generated with optool) is in ['models/dust.py'](./models/dust.py).

Specific functions to write radmc files are in ['models/write_radmc_files.py'](./models/write_radmc_files.py).


### Plotting
Some convenience plotting functions can be found in ['models/make_plots.py'](./models/make_plots.py).

### Chemical Model
The chemical model files are written based off of the radmc outputs and re-gridded for the chemical code. The setups for the chemistry run are based off of template values and parameters set during the model setup. The details of how this is performed are in ['models/prepchem/chem_io.py'](./models/prepchem/chem_io.py').
