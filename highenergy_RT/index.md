### High Energy Radiative Transfer
[back to main](../README.md)

This is a list of the radiative processes done before the chemistry. 

## UV contribution

Original Bethell code based on [Bethell and Bergin 2011](https://iopscience.iop.org/article/10.1088/0004-637X/739/2/78/meta#apj400374s3)

Two components (from original paper)

1. fUV continuum from stellar source with scattering computed using Henyey–Greenstein phase function. (this is the main MC done in the code)
2. Lyman alpha transport (done afterwards based on template spectrum?)

## X-ray contribution

X-ray Bethell code [Bethell and Bergin 2011b](https://iopscience.iop.org/article/10.1088/0004-637X/740/1/7/meta)

X-ray fluxes are computed similarly, with some template set of energies and x-ray opacities (**missing this file**) and output into ```.dat``` files

Additional component (unclear where this goes): ```ZetaHe``` and ```ZetaH2``` - finds ionization rates and outputs them. Is this used in the chem?

**Note**

In principle, both of these could be done with the ```radmc3d mcmono ``` method either by adopting simple continuum models or using template spectra for the fUV and X-ray for which radmc3d would provide mean intensity fields for the UV and X-ray at every point.

UV Henyey-Greenstein for ```radmc-3d mcmono``` should be done with ```scattering_mode=2```

Then the gas temperature could be recalculated just like in ```torus2chem.py```
This has the advantage that the model wouldn't be downsampled and interpolated for the high energy RT, so would be consistent with the values put into the original radmc model. (really cuts down on the amount of interpolation in general if this is computed on the original spherical grid)

+[radmc3d mcmono](https://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/manual_radmc3d/dustradtrans.html#sec-dust-monochromatic-monte-carlo)

The final values for the temperature could be computed on the original grid and then interpolated with everything else to the r,z slices for the chemistry.

Also has the advantage of not having to have IDL to make this work...

See Appendix in this [DIANA paper](https://www.aanda.org/articles/aa/pdf/2016/02/aa26538-15.pdf) for fUV and X-ray continuum models. 

## Cosmic ray flux 
Attenuation in a column down the midplane.

## Interstellar Radiation Field (ISRF)
Radmc-3D can take an ISRF to add to the dust temperature computation. Since the chemistry takes it as well, it makes sense to be consistent with the dust temperature calculations and use the same one for each. 


