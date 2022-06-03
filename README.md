# wedgeford
Code repository for the Warm Emission Due to Gas Envelopes Falling Onto Rotating Disks (WEDGEFORD) project


## Basic Code Workflow 

+ Create a physical model of a disk and infalling envelope for some number of dust species:
    + Grid
        + spherical coordinate grid in `r,theta,phi`
            + `model.make_grid()` (3d version)
            + `model.r, model.theta, model.phi` (1d versions)
            + output file: `"amr_grid.inp"`  
    + Gas component
        + Ulrich-envelope + gas disk with thermal disk scale height set by irradiation
        + output: array of 3d volume densities 
            + `model.rho_embedded(fluid=0)`  
    + Dust
        + Several dust species, with sizes and opacities generated from library
        + output: array of 3d volume densities
            + `model.rho_embedded(fluid=1,2...etc)`    
+ Run `radmc3d mctherm` on model outputs
    + output file: `"dust_temperature.dat"` 
    + output: 3-d dust temperature on the original spherical grid 
+ Slice a quadrant in r-z and interpolate (densities and temperatures) onto a cylindrical grid in r and z
    + `model.make_rz_uniform()`
    + `model.make_quadrant(,fluid=...)`
    + output file: `gas..out`
+ Compute cosmic ray attentuation based on column density at each point in the r,z grid
    + `model.zeta()`
+  Run UV and X-ray radiative transfer
    + input: `gas..out` 
    + output: UV and X-ray fluxes at each point in the r,z grid
    + output files: `uv..dat` and `xray.dat`
+ Generate input files for chemical models  
    + chemistry inputs: `uv..dat`, `xray..dat`, `1_env..inp` 
    + `taurus_2_chem.py` makes `1_env..inp`
+ Run chemical models  


## Code components
### [Tutorial notebook](./make_model.ipynb): basic notebook to set up a problem in radmc3d.

+ [models](./models/index.md)
  + Creation of model object
  + Radiative Transfer of model to calculate dust temperature
  + Cosmic ray attenuation
+ [chemistry](./chemistry/index.md)
  + Full chemical modeling

