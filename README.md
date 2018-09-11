
# Finding Emission Lines In Spectra (FELIS)

README for the template matching software FELIS.

If you find FELIS useful please reference this FELIS GitHub repository at https://github.com/kasperschmidt/FELIS. 

## Table of Content
<a href="FELISlogo.png"><img src="FELISlogo.png" align="right" height="180" ></a>

- [Description](#description)
- [Installing FELIS](#installing-felis)
- [Script Overview](#script-overview)
- [Running FELIS](#running-felis)
- [FELIS Output](#felis-output)

## Description

The small tool for Finding Emission Lines In Spectra (FELIS) is build in Python, and can be used to search for spetral features like emission lines in extracted 1D spectra.

## Installing FELIS

FELIS does not require any installation. Simply cloning the FELIS GitHub repository and importing the scripts should be enough.
Hence, FELIS is "installed" by doing:
```
git clone https://github.com/kasperschmidt/FELIS.git
```
After adding the FELIS directory to the `PYTHONPATH` or changing location to the FELIS directory, FELIS can be imported in `python` with:
```python
import felis
```
If the import does not generate any errors, i.e., if the FELIS dependencies and requirements are met, FELIS is ready for use. Apart from `astropy` FELIS only uses standard Python packages that can be installed with `pip`.


## Script Overview

The template matching itself is handled by `felis.py`, whereas `felis_build_template.py` an be used to generate the tempaltes to match to the 1D spectra. Below, a few of the functions and tourines in these two scripts are described, but for futer details, refer to the header of each routine.

- `felis.match_templates2specs()`
  - The wrapper that carries out the matching of tempaltes with spectra, and returns a dictionary with the results, including the estimated S/N of each template match.
  
- `felis.cross_correlate_template()`
  - The function that performs the cross-correlation, or rather the template match, between each tempalte and spectrum provided, by minimizig the chi2 between the template and the input data.

- `felis_build_template.build_template()`
  - The main function to generate templates that can be matched to the 1D spectra with .

## Running FELIS

In the following a few useful examples of how to produce outputs with FELIS are given. For examples on how to run individual pieces of code, please refer to the individual code headers.

A standard run of felis, matching a collection of templats to a single spectrum (which is first put on the FELIS/TDOSE format) can be done with the following line of commands:
```python
import felis, glob

inputspec  = pyfits.open(spec)[1].data
wave    = inputspec['Wavelength']
flux    = inputspec['Spectrum']
fluxerr = inputspec['Error']

felisspec   = /path/to/location/of/FELISformatspec.fits
templatedir = /dircotiry/containing/templates/to/match/to/

print(' - Saving spectrum in the FELIS/TDOSE format: \n   '+outfile)
    felis.save_spectrum(felisspec,wave,flux,fluxerr,headerinfo=None,overwrite=True,verbose=verbose)

print(' - Grabbing templates from '+templatedir)
temps   = glob.glob(templatedir+'basename_of_templates*.fits'')

plotdir  = '/directory/where/plots/are/stroed/to/'
print(' - Matching templates to spectrum using FELIS')
picklefile = '/location/of/output/picklefile/felisresults.pkl'

ccdic = felis.match_templates2specs(temps, [felisspec],[2.35], picklefile, wavewindow=[15], plotdir=plotdir, wavecen_restframe=[5007.0], vshift=[0], min_template_level=1e-4, plot_allCCresults=True, subtract_spec_median=True)
```

Templates can be generated with expression like the following:
```python
import felis_build_template as fbt
import numpy as np

tcdic = {}
tcdic['D1900']                  = ['DELTA', 1900.0, 10.0,           'Delta function at 1900A']
tcdic['CIII1']                  = ['GAUSS', 1907.0, 0.5, 0.0, 10.0, 'CIII]1907A']
tcdic['CIII2']                  = ['GAUSS', 1909.0, 0.5, 0.5, 5.0,  'CIII]1909A']
tcdic['CONT']                   = ['CONT', 1.0, -0.03, 1908.0,       'Continuum with flux 1.0 at 1908 + slope -0.03']
break_wave                      = np.arange(1800,2000,2.0)
break_flux                      = break_wave*0.0
break_flux[break_wave > 1920.0] = 5.0
tcdic['B1920']                  = ['FEATURE',break_wave, break_flux,'Break at 1920 Angstrom']

fbt.build_template([1870,1980,0.1],tcdic,tempfile='./felis_template_example.fits',overwrite=True)
```
The tempalte generated with the above command is stored in `felis_template_example.fits`.

## FELIS Output

A template match carried out with FELIS generates a dictionary containng containing all the cross-correlation results. The dictionary will be returned directly but also saved to disk as a pickled filed. This file can be loaded with: `felis.load_picklefile(picklefilename)`

The returned dictionary keys contain the individual spectra with corresponding template match results. Hence, `dictionary.keys()` is the list of spectra that have been crossmatched.

Each entry in the dictionary (`dictionary[key]`) contains the following entries:

 - `wavelengths`: The wavelength vector used for the cross-correlation matching of each of the N templates
 - `templatevec`: List of each of the N templates matched to spectrum
 - `zspec`: Spectroscopic redshift for spectrum if known.
 - `zCCmaxvec`: The redshift corresponding to max S/N for each of the N templates matched.
 - `ccresultsarray_flux`: Flux vectors for each of the N templates matched
 - `ccresultsarray_variance`: Variance vecotrs for each of the N templates matched
 - `ccresultsarr_S2N`: S/N vector for each of the N templates matched
 - `ccresultsarr_chi2`: Chi squared values for the cross-correlation of each of the N templates matched
 - `ccresultsarr_Ngoodent`: The number of good pixels used in the cross correlation for each of the N templates matched
 - `S2NCCmaxvec`: Vector with max(S/N) values for N templates matched
 - `continuumlevel`: The 'continuum level' of the spectrum removed in the cross-correlation. Currently the value is simply the median of the spectrum.
 - `vshift`: If a velocity shift was provided for the template match this is stored here

The results picklefile can be used to assemble sub-sample results (e.g., S/N cuts) based on the template cross-correlations with `felis.selection_from_picklefile()` and individual entries can be plotted using `felis.plot_picklefilecontent()`.