# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
import felis
import felis_build_template as fbt
import matplotlib.pyplot as plt
from astropy.io import fits
import pdb
import sys
import scipy
import numpy as np
import collections

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def build_template(wavelengths,templatecomponents,noise=None,
                   tempfile='./felis_template_RENAME_.fits',
                   plottemplate=True,zoomxplot=None,overwrite=False,verbose=True):
    """
    Build a spectral template based on a dictionatry specifying the tempalte components
    like emmision line and continuum shape.

    --- INPUT ---
    wavelengths                 Wavelength range and resolution [lambda_min, lambda_max, delta_lambda] to build
                                template over. Any template component centered outsize this range is ignored.
    templatecomponents          A dictionary specifying the components of the template to build
                                The inputs of the dictionary could be:

                                An emission line:
                                    To add components to the spectrum chose one of the following setups, using the
                                    first entry of the keyword list to indicate the component type such that:
                                        tcdic['keyname'] = [comptype, param1, param2, ...]
                                        NB: Keep 'keyname' to <=5 characters to have a maximum fits header
                                            lkeyword ength of 8 characters (3 are added when saving the info)
                                    Available component types are:
                                    ['GAUSS', mean, sigma, skew, scaling, 'info'] # Gaussian (Skewed) profile
                                    ['DELTA', position, , total flux, 'info']     # Delta function. Will be added at
                                                                                  # wavelength nearest postion in the
                                                                                  # generated template
                                    ['LSF', sigma, 'info']                        # Gaussian LSF to convolve template with
                                    ['CONT', level, slope, lam0, 'info']          # linear continuum on the form
                                                                                  # C = level + slope * (lam-lam0)
                                                                                  # Alternatively, a pre-defined continuum
                                                                                  # Can be aded as a 'flux feature'
                                    ['FEATURE',wavelength,flux,'info']            # A spectral feature (like breaks and
                                                                                  # broad odd line pofiles) defined by a
                                                                                  # wavelength and flux vector

                                In the above, 'info', is a string with info to write to template fits header
                                when listing the template components
    noise                       To add noise to the template provide one of the following:
                                    ['POISSON',mean]       Drawing noise from Poisson distribution around mean
                                    ['GAUSS',mean,sigma]   Drawing noise from Gaussian distributed error spectrum
                                    ['CONSTANT',value,]    Drawing noise from constant error spectrum
                                    ['SPECTRUM',wave,flux] Error spectrum to draw noise from after interpolation                                                                      to template range
    tempfile                    Name of fits file to store final template to
    plottemplate                Plot the generated template?
    zoomxplot                   Only show plot in zoom region given as [lambda_min,lambda_max]
    overwrite                   Overwrite existing template?

    --- EXAMPLE OF USE ---
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

    tcdic = {}
    tcdic['CONT']                   = ['CONT',  1.0, 0.0, 0.0,          'Flat continuum at 1.0']
    tcdic['CIII1']                  = ['GAUSS', 1907.0, 0.5, 0.0, 10.0, '[CIII]1907A']
    tcdic['CIII2']                  = ['GAUSS', 1909.0, 0.5, 0.0, 5.0,  'CIII]1909A']
    noise                           = ['GAUSS', 2.0, 1.0]

    fbt.build_template([1900,1920,0.1],tcdic,tempfile='./felis_template_CIIIdoublet.fits',noise=noise,overwrite=True)


    """
    if verbose: print(' - Setting up template output ')
    headerdic = collections.OrderedDict()
    wavevec   = np.arange(wavelengths[0],wavelengths[1],wavelengths[2])
    fluxvec   = wavevec*0.0
    LSF       = False  # default is no LSF convolution
    N_LSF     = 0      # count number of LSF setups provided - will only handle 1 set of parameters.

    for key in templatecomponents.keys():
        if templatecomponents[key][0].upper() == 'GAUSS':
            param     = templatecomponents[key][1], templatecomponents[key][2], templatecomponents[key][3]
            gaussvec  = fbt.gauss_skew(wavevec,*param)
            fluxvec   = fluxvec + templatecomponents[key][4] * gaussvec
            for kk in [0,1,2,3,4]:
                headerdic['F'+key+'_'+str(kk)] = [templatecomponents[key][kk],templatecomponents[key][-1]]

        elif templatecomponents[key][0].upper() == 'CONT':
            fluxvec  = fluxvec + templatecomponents[key][1] + \
                       templatecomponents[key][2]*(wavevec - templatecomponents[key][3])
            for kk in [0,1,2]:
                headerdic['F'+key+'_'+str(kk)] = [templatecomponents[key][kk],templatecomponents[key][-1]]

        elif templatecomponents[key][0].upper() == 'DELTA':
            dwave    = np.abs(wavevec - templatecomponents[key][1])
            waveent  = np.where(dwave == np.min(dwave))[0]
            fluxvec[waveent] = fluxvec[waveent] + templatecomponents[key][2]
            for kk in [0,1,2]:
                headerdic['F'+key+'_'+str(kk)] = [templatecomponents[key][kk],templatecomponents[key][-1]]

        elif templatecomponents[key][0].upper() == 'LSF':
            key_LSF = key
            LSF     = True
            N_LSF   = N_LSF + 1

        elif templatecomponents[key][0].upper() == 'FEATURE':
            func       = scipy.interpolate.interp1d(templatecomponents[key][1],templatecomponents[key][2],
                                                    kind='linear',fill_value=0.0)
            flux_feat = func(wavevec)
            fluxvec    = fluxvec + flux_feat
            headerdic['F'+key+'_1'] = [templatecomponents[key][0],templatecomponents[key][-1]]

        elif templatecomponents[key][0].upper() == 'SKIP':
            continue

        else:
            sys.exit('Invalid template component "'+templatecomponents[key][0]+'"')

    if noise is not None:
        if noise[0] is 'POISSON':
            if verbose: print(' - Drawing Poissonian noise')
            fluxerr  = np.abs(np.random.poisson(noise[1],fluxvec.shape).astype(float))
            headerdic['FNOISE_1'] = [noise[1],'Mean of Poissonian noise']
        elif noise[0] is 'CONSTANT':
            if verbose: print(' - Drawing from constant noise level')
            fluxerr = fluxvec*0.0 + noise[0]
            noisepix  = fluxerr*0.0
            for nn, noisesigma in enumerate(fluxerr):
                noisepix[nn] = np.random.normal(0.0, noisesigma, 1)

        elif noise[0] is 'GAUSS':
            if verbose: print(' - Adding Gaussian noise')
            headerdic['FNOISE_1'] = [noise[1],'Mean of Gaussian noise']
            headerdic['FNOISE_2'] = [noise[2],'Sigma of Gaussian noise']
            fluxerr = np.abs(np.random.normal(noise[1], noise[2], fluxvec.shape))
            noisepix  = fluxerr*0.0
            for nn, noisesigma in enumerate(fluxerr):
                noisepix[nn] = np.random.normal(0.0, noisesigma, 1)

        elif noise[0] is 'SPECTRUM':
            if verbose: print(' - Adding noise based on interpolation of provided noise spectrum ')
            noisewave = noise[1]
            noiseflux = noise[2]
            noisefunc = scipy.interpolate.interp1d(noisewave,noiseflux,kind='linear',fill_value=0.0)
            fluxerr   = noisefunc(wavevec)

            noisepix  = fluxerr*0.0
            for nn, noisesigma in enumerate(fluxerr):
                noisepix[nn] = np.random.normal(0.0, noisesigma, 1)

            headerdic['FNOISE_1'] = [np.mean(noiseflux),'Average noise of provided noise spectrum ']
        else:
            sys.exit('Invalid noise component "'+noise[0]+'"')
    else:
        fluxerr  = fluxvec*0.0
        noisepix = fluxerr
        headerdic['FNOISE_1'] = ['None','No noise added']

    fluxvec = fluxvec + noisepix

    headerdic['FERR_1'] = ['noise','Uncertainty on flux = noise from FNOISE keys']

    if verbose: print(' - Checking for LSF parameters ')
    if LSF:
        if N_LSF == 1:
            if verbose: print('   Convolving template with LSF setup found in component dictionary: ')
            LSFparam = templatecomponents[key_LSF]  #['LSF', sigma, 'info'] # Gaussian LSF to convolve template with
            if verbose: print('   LSF info: '+LSFparam[2])
            if verbose: print('   sigma_LSF = '+str(LSFparam[1]))

            dwave       = np.median(np.diff(wavevec))
            Nsigma      = 5.0
            windowwidth = np.int( (LSFparam[1]/dwave) * 2.0 * Nsigma ) # window width of +/-10 sigma

            if windowwidth % 2 == 0: # window width is even.
                windowwidth = windowwidth + 1

            if windowwidth < 5:
                windowwidth = 5
                print('   LSF sigma resulted in a LSF convolution window of less than 5 wavelength pixels. '
                      'Using 5 pixels for LSF gaussian.')

            if windowwidth > len(wavevec):
                sys.exit('LSF sigma +/- '+str(Nsigma)+'sigma restults in a LSF Gaussian '
                         'convolution kernel larger than considered wavelength range ('+str(windowwidth)+
                         ' > '+str(len(wavevec))+')')

            wave_LSF    = wavevec[0:windowwidth]
            param       = np.mean(wave_LSF[int(windowwidth/2.)]), LSFparam[1], 0.0
            gauss_LSF   = fbt.gauss_skew(wave_LSF,*param) * dwave # normalize Gaussian

            # Check that gaussian is normalized ...
            # print(' ----> Trapetzoidal integration of LSF guass = '+str(scipy.integrate.trapz(gauss_LSF,wave_LSF)))
            # print(' ----> sum of LSF guass (PDF sums to 1)      = '+str(scipy.integrate.trapz(gauss_LSF,wave_LSF)))
            # plt.plot(wave_LSF,gauss_LSF,'bo')
            # plt.savefig('./LSFgenerated.pdf')
            # plt.clf()

            fluxvec   = np.convolve(fluxvec, gauss_LSF, 'same')#
        else:
            sys.exit(' Found '+str(N_LSF)+' LSF setups when building FELIS template; only one LSF setup can be provided')
    else:
        if verbose: print('   None found - no convolution with LSF ')

    if verbose:
        print('   - Spectrum stats:    flux(std,mean,median,max,min)'+
              str([np.std(fluxvec),np.mean(fluxvec),np.median(fluxvec),np.min(fluxvec),np.max(fluxvec)]))
        print('                     fluxerr(std,mean,median,max,min)'+
              str([np.std(fluxerr),np.mean(fluxerr),np.median(fluxerr),np.min(fluxerr),np.max(fluxerr)]))
        print('                         S2N(std,mean,median,max,min)'+
              str([np.std(fluxvec/fluxerr),np.mean(fluxvec/fluxerr),np.median(fluxvec/fluxerr),np.min(fluxvec/fluxerr),np.max(fluxvec/fluxerr)]))

    felis.save_spectrum(tempfile,wavevec,fluxvec,fluxerr,headerinfo=headerdic,overwrite=overwrite,verbose=verbose)

    if plottemplate:
        fbt.plot_template(tempfile,showerr=True,zoomx=zoomxplot,verbose=verbose)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def gauss_skew(x, *p):
    """
    Return the PDF of a skewed normal distribution where alpha is the skewness paramter
    alpha = 0 represent a non-skewed regular Gaussian with mean mu and standard deviation sigma

    --- INPUT ---
    x    Wavelength vector to calculate gaussian for
    *p   The parmaters of the skewed Gaussian PDF, i.e., p = mu, sigma, alpha

    """
    mu, sigma, alpha = p

    t = (x-mu)/sigma
    A = 1.0 / ( sigma*np.sqrt(2.0*np.pi) )

    gauss_pdf = A * np.exp(-t**2.0 / 2.0)
    gauss_cdf = (1.0 + scipy.special.erf(alpha*t/np.sqrt(2.0))) / 2.

    gaussskew_pdf = 2.0 * gauss_pdf * gauss_cdf

    return gaussskew_pdf

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_template(templatefits,zoomx=None,showerr=True,verbose=True):
    """
    Plotting template

    --- INPUT ---
    templatefits     Fits file containg binary table with template
    showerr          Show the errorbars on the flux values
    verbose          Toggle verbosity

    --- EXAMPLE OF USE ---
    templatefits = './felis_template_TEST.fits'
    fbt.plot_template(templatefits,showerr=True,verbose=True)

    """
    specdata = fits.open(templatefits)[1].data
    datahdr  = fits.open(templatefits)[1].header
    plotname = templatefits.replace('.fits','_plot.pdf')
    if verbose: print(' - Setting up and generating plot')
    fig = plt.figure(figsize=(7, 5))
    fig.subplots_adjust(wspace=0.1, hspace=0.3,left=0.1, right=0.97, bottom=0.10, top=0.95)
    Fsize    = 12
    lthick   = 2
    marksize = 4
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=Fsize)
    plt.rc('xtick', labelsize=Fsize)
    plt.rc('ytick', labelsize=Fsize)
    plt.clf()
    plt.ioff()

    try:
        plt.xlabel(' Wavelength ['+datahdr['TUNIT1']+']')
    except:
        plt.xlabel(' Wavelength [?]')

    xvalues  = specdata['wave']
    yvalues  = specdata['flux']
    yvalues2 = specdata['s2n']
    yerr     = specdata['fluxerror']

    if zoomx is not None:
        goodent  = np.where( (xvalues > zoomx[0]) & (xvalues < zoomx[1]) )
        if len(goodent) == 0:
            sys.exit('No data in zoomx region to plot')
        xvalues  = xvalues[goodent]
        yvalues  = yvalues[goodent]
        yvalues2 = yvalues2[goodent]
        yerr     = yerr[goodent]

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plt.subplot(2,1,1)
    plt.plot(xvalues,yvalues,'go',lw=lthick, markersize=marksize,alpha=1.0,)
    if showerr:
        plt.fill_between(xvalues, yvalues-yerr, yvalues+yerr,
                         color='green',alpha=0.5)

    try:
        plt.xlabel(' Wavelength ['+datahdr['TUNIT1']+']')
    except:
        plt.xlabel(' Wavelength [?]')

    try:
        plt.ylabel(' Flux ['+datahdr['TUNIT2']+']')
    except:
        plt.ylabel(' Flux [?]')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plt.subplot(2,1,2)
    plt.plot(xvalues,yvalues2,'go',lw=lthick, markersize=marksize,alpha=1.0,)
    try:
        plt.xlabel(' Wavelength ['+datahdr['TUNIT1']+']')
    except:
        plt.xlabel(' Wavelength [?]')
    plt.ylabel(' S/N ')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print('   Saving plot to '+plotname)
    plt.savefig(plotname)
    plt.clf()
    plt.close('all')
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
