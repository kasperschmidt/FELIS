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
