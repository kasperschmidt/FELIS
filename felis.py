# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
import felis
from astropy.io import fits
import pdb
import time
import os
import sys
import glob
import numpy as np
import collections
import astropy
import collections
from astropy import wcs
from astropy.coordinates import SkyCoord
from astropy import units
import astropy.convolution
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def find_emission_lines_in_spectra(linelist,continuum):
    """
    Wrapper for felis.find_emission_lines_in_spectrum

    --- INPUT ---

    --- EXAMPLE OF USE ---
    import felis
    spectra = glob.glob('./tdose_spectrum*.fits')
    felis.build_template(spectra,verbose=True)

    """

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def find_emission_lines_in_spectrum(spec_wace_spec_flux,template_library,verbose=True):
    """
    Finding emission lines in spectrum, by cross-correlating a template library with the input spectrum

    --- INPUT ---

    --- EXAMPLE OF USE ---
    import felis

    """
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def build_template_library(linelist,continuum):
    """
    Assemble a template library for cross correlation using templates buld with felis.build_template()

    --- INPUT ---

    --- EXAMPLE OF USE ---
    import felis

    """
    templib = {}

    return templib
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def build_template(linelist,continuum):
    """
    Build a spectral template based on an input line list and continuum

    --- INPUT ---

    --- EXAMPLE OF USE ---
    import felis

    felis.build_template(verbose=True)

    """
    if verbose: print(' - Building template')


    if verbose: print(' - Storing template')


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def cross_correlate(spectrum,template,continuum):
    """
    Build a spectral template based on an input line list and continuum

    --- INPUT ---

    --- EXAMPLE OF USE ---
    import felis

    """

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def build_mock_spectrum(linelist,lineratios,wavelength,continuum,noise='Gauss',verbose=True):
    """
    Function to build a mock spectrum with know lines, flux ratios, continuum and noise properties.

    --- INPUT ---

    --- EXAMPLE OF USE ---
    import felis

    """

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def save_spectrum(outfile,wave,flux,fluxerr,headerinfo=None,overwrite=False,verbose=True):
    """
    Function saving a spectrum (or template) to a fits file

    --- INPUT ---
    outfile         Name of output file to generate
    wave            Wavelength in units of Angstrom
    flux            Flux to store to fits file
    fluxerr         Error un flux
    headerinfo      To add info to the header provide it to this keyword as a dictionary on the format:
                       headerinfo = {'KEYNAME1':[VALUE1,INFOCOMMENT1], 'KEYNAME2':[VALUE2,INFOCOMMENT2]}
    overwrite       Overwrite existing file?
    verbose         Toggle verbosity

    --- EXAMPLE OF USE ---
    import felis
    import numpy as np
    wave        = np.arange(1,100,1)
    flux        = wave**0.5 + (wave*0.0+2.3)
    fluxerr     = np.sqrt(flux)
    headerinfo  = {'FEL_EL1':[1548,'Wave of FELIS line 1 in template'], 'FEL_EL2':[1551,'Wave of FELIS line 2 in template']}
    outfile     = './spectrum_output.fits'

    felis.save_spectrum(outfile,wave,flux,fluxerr,headerinfo=headerinfo,overwrite=False,verbose=True)


    """
    S2N = flux/fluxerr

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Saving wavelenght and flux values to \n   '+outfile)
    mainHDU = fits.PrimaryHDU()       # primary HDU
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    c1 = fits.Column(name='wave',      format='D', unit='ANGSTROMS', array=wave)
    c2 = fits.Column(name='flux',      format='D', unit='', array=flux)
    c3 = fits.Column(name='fluxerror', format='D', unit='', array=fluxerr)
    c4 = fits.Column(name='s2n',       format='D', unit='', array=S2N)

    coldefs = fits.ColDefs([c1,c2,c3,c4])
    tbHDU   = fits.BinTableHDU.from_columns(coldefs) # creating default header

    # writing hdrkeys:'---KEY--',                             '----------------MAX LENGTH COMMENT-------------'
    tbHDU.header.append(('EXTNAME ','SPEC1D'                     ,'cube containing source'),end=True)
    if headerinfo is not None:
        for key in headerinfo.keys():
            tbHDU.header.append((key,headerinfo[key][0],headerinfo[key][1]),end=True)

    hdulist = fits.HDUList([mainHDU,tbHDU])
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    hdulist.writeto(outfile, overwrite=overwrite)
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
