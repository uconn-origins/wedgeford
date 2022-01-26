### High Energy Radiative Transfer
[back to main](../README.md)

This is a list of the radiative processes done before the chemistry. 

## UV contribution

Original Bethell code based on [Bethell and Bergin 2011](https://iopscience.iop.org/article/10.1088/0004-637X/739/2/78/meta#apj400374s3)

Two components (from original paper)

1. fUV continuum from stellar source with scattering computed using Henyey–Greenstein phase function. (this is the main MC done in the code)
2. Lyman alpha transport --> from original Fogel Paper: compute the H2 fraction from ionization rate
![image](https://user-images.githubusercontent.com/20684970/151234729-5cc5d3a3-f13e-4b4c-9c8e-43207e4f8bab.png)
then recompute isotropic scattering from Lyman alpha photons with atomic hydrogren cross-section

*Note on Lyman alpha (from reading too much about this, god help me)*
the Bethell and Bergin MCMC computes the location of the H2/H dissociation layer iteratively (gross). Really only the bottom surface is really important because it's the scattering surface that scatters Lya photons deeper into the disk, thus increasing UV photon depth. The bottom surface is pretty much the boundary between H2/H, which is sort of determined by the amount of self-shielding you have. basically at a certain column, you have enough self-shielding that H2 doesn't dissociate. 

Anyway, it seems to me that the way around the iterative H2/H nonsense, is just to find the column where self-shielding matters, define that as the scattering surface and set up radiative transfer with a population of fake "dust" H atoms (at 0.5 of the gas density) that have the scattering opacities of the lyman alpha line you're putting in to the code and use the isotropic scattering mode (will be very fast) and see how far the UV penetrates.

Benefits: This will pretty much give you the maximum UV field penetration depth without having to worry about iterating over where your dust at etc because its built in! Aaaaaand has the added benefit that we can compare to the fUV + Lyman line without extra resonant scattering and see if it even matters for this thing, like an upper/lower limit kinda deal. 

I prefer this phenomenology because as B+B 2011 admits, their iterative photodissociation stuff isn't even self-consistent/fully correct anyway and we can quickly see if this applies to the problem at all. Because if the scattering layer occurs deep beneath where UV can't even go to begin with ( which may be the case here) then we don't even gotta worry about it at all. (Also other people just kinda admit on whiffing on the whole finding H layer to begin with and make even more egregious estimates, so I think we're fine). 

Second, I just read way too many CTTS spectroscopy papers just to learn that no one (except maybe like only a wee bit Nuria) even understands the physical production of Lya in these systems other than "waves hands accretion, waves hands outflows, we're all just waving our hands" -- like I don't think people even know if this is applicable to Class 0/I sources to begin with (LoooooooooL). Maybe ~*~veiling~*~ saves the day (also lol).

So I think maybe running with 2 edge cases is the smartest way around this. 


## X-ray contribution

X-ray Bethell cross-sections [Bethell and Bergin 2011b](https://iopscience.iop.org/article/10.1088/0004-637X/740/1/7/meta)

X-ray fluxes are computed similarly, with some template set of energies and x-ray opacities (**missing this file**) and output into ```.dat``` files

X-ray cross-sections are computed with fitting functions found in original table in Bethell and Bergin (given in per H nucleus)

![image](https://user-images.githubusercontent.com/20684970/150820632-184e346a-96db-498f-9415-0bfd225391d4.png)

and components are combined:

![image](https://user-images.githubusercontent.com/20684970/150821591-ef28bf2c-47d5-4a73-9a38-13167814d1f5.png)

with some assumption of the "blanketing factor" and "dust settling" - both of these just scale the contribution to the cross-section from dust.
The xray-radiative transfer in the Bethell code uses these for the local MC calculations inputting an epsilon for each point in the local grid.

The ionization rates for He and H are computed based on the X-ray flux at each location. 


**Note**

In principle, both of these could be done with the ```radmc3d mcmono ``` method either by adopting simple continuum models or using template spectra for the fUV and X-ray for which radmc3d would provide mean intensity fields for the UV and X-ray at every point.

UV Henyey-Greenstein for ```radmc-3d mcmono``` should be done with ```scattering_mode=2```
X-ray ```radmc-3d mcmono``` can be done with the absorption coefficients scaled. (i.e. for some constant/average dust to gas ratio of the small dust to include the constant gas component there as well). We don't really need to do the RT with the large grains in this case as the absorption coefficients are much much smaller for larger grains and will affect the opacities a lot less. 

These opacities can be added into the ```dust_kappa...inp``` files based on the chosen model. (```mc therm``` will run fine as the wavelengths will be limited to the stellar spectrum and ISRF). 

Then the gas temperature could be recalculated just like in ```torus2chem.py``` on the original grid. 

This has the advantage that the model wouldn't be downsampled and interpolated for the high energy RT, so would be consistent with the values put into the original radmc model. (really cuts down on the amount of interpolation in general if this is computed on the original spherical grid)

+[radmc3d mcmono](https://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/manual_radmc3d/dustradtrans.html#sec-dust-monochromatic-monte-carlo)

The final values for the temperature could be computed on the original grid and then interpolated with everything else to the r,z slices for the chemistry.

Also has the advantage of not having to have IDL to make this work...

See Appendix in this [DIANA paper](https://www.aanda.org/articles/aa/pdf/2016/02/aa26538-15.pdf) for some simple fUV and X-ray continuum models.

Then can use the same output writer in the compilation scripts to generate the ```.dat``` files for the chemistry.

## Cosmic ray flux 
Attenuation in a column down the midplane based on model fits from Cleeves et al 2013.

## Interstellar Radiation Field (ISRF)
Radmc-3D can take an ISRF to add to the dust temperature computation. Since the chemistry takes it as well, it makes sense to be consistent with the dust temperature calculations and use the same one for each. 



