# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
import felis
from astropy.io import fits
import scipy.signal
import pdb
import matplotlib.pyplot as plt
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
def find_emission_lines_in_spectrum(spec_wace,spec_flux,template_library,verbose=True):
    """
    Finding emission lines in spectrum, by cross-correlating a template library with the input spectrum

    --- INPUT ---

    --- EXAMPLE OF USE ---
    import felis

    """
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def get_template_library(linelist,continuum):
    """
    Assemble a library of templates (build with felis_build_template.py) to crosscorrelate spectrum with

    --- INPUT ---

    --- EXAMPLE OF USE ---
    import felis

    """
    templib = {}

    return templib

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def minimize_chi2(data,data_sigma,template,verbose=True):
    """
    Minimizing chi2 for a 1D tempalte matching and returning the corresponding scaling, alpha.

    Done using that
        chi**2         = Sum_i (d_i - alpha * t_i)**2 / sigma_i**2
    and that dchi**2/dalpha = 0 implies
        alpha          = ( Sum_i d_i*t_i/sigma_i**2 ) / ( Sum_i t_i**2 / sigma_i**2 )
        sigma_alpha**2 = 1 / ( Sum_i t_i**2 / sigma_i**2 ) = alpha_variance

    Here, d is the data (flux from spectrum), sigma is the uncertainy on the data, t is the template,
    alpha is the flux scaling of the template, sigma_alpha is the uncertainty on alpha, and i runs
    over the pixels in the spectrum.

    --- TEMPALTE ---

    data            Data to match template to
    data_sigma      Uncertainty on data, i.e., sqrt(variance)
    template        Template to serach for in data. Should be of the same lenght as data
    verbose         Toggle verbosity

    """
    if len(data) == len(template):
        Npix = len(data)
    else:
        sys.exit('The length of the data and template should be the same; it is not.')
    if verbose: print(' - Will calculate and minimize chi**2 for data and template of length '+str(Npix))
    goodent  = np.where((data_sigma > 0) & (np.isfinite(data)) & (template != 0))[0]
    #goodent  = np.where((data_sigma > 0) & (data > 0) & (template != 0))[0]
    Ngoodent = len(goodent)

    if Ngoodent == 0:
        if verbose: print(' - No entries left where data_sigma > 0 & data is finite & tempalte != 0, so returning 0s')
        chi2_min        = 0.0
        alpha           = 0.0
        alpha_variance  = 0.0
        S2N             = 0.0
    else:
        alpha          = np.sum( data[goodent] * template[goodent] / data_sigma[goodent]**2 ) / \
                         np.sum( template[goodent]**2 / data_sigma[goodent]**2 )

        alpha_variance = 1.0 / np.sum( template[goodent]**2 / data_sigma[goodent]**2 )
        S2N            = alpha / np.sqrt(alpha_variance)
        chi2_min       = np.sum( (data[goodent] - alpha * template[goodent])**2 / data_sigma[goodent]**2 )

    return alpha, alpha_variance, S2N, chi2_min, Ngoodent
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def cross_correlate_template(spectrum,template,z_restframe=None,waverange=None,plotresult=None,
                             min_template_level=1e-2,spec_median_sub=False,verbose=True):
    """


    --- INPUT ---
    spectrum            fits spectrum to find a (cross-correlation) match to template for
    template            fits template to correlate with
    z_restframe         To perform cross-correlation in rest-frame (shifting the spectrum) provide a redshift.
    waverange           To prevent running the cross-correlation on the full spectrum, provide the wavelength
                        range (based on wavelength-vector in spectrum) to restrict the correlation to.
    plotresult          Provide an output file name to illustrate the location for the best
                        correlation between the spectrum and template
    min_template_level  The template is interpolated to the wavelength grid of the spectrum and extrapolated
                        beyond it's edges if nescessary. In this extrapolation (assuming the template goes to ~0
                        at the edges), very small values (e.g., <1e-20) can be returned. To set these to 0.0
                        provide a level below which all values in the interpolated template are treated as 0s.
    spec_median_sub     If True, the median level (continuum?) of the spectrum is subtracted prior to matching.
    verbose             Toggle verbosity

    --- EXAMPLE OF USE ---
    import felis

    specdir  = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/felis_testing/'
    spectrum = specdir+'uves_felis_mock_MUSEspectrum_noisesigma1p0.fits'
    spectrum = specdir+'uves_felis_mock_MUSEspectrum_noisesigma0p05.fits'
    template = specdir+'uves_felis_template_CIIIdoublet_sig_0p25_fluxCIII1_2p0_fluxratio_0p5.fits'
    plotname = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/felis_testing/test_CCResultPlot.pdf'

    s_wave, ccresults, max_S2N, max_z = felis.cross_correlate_template(spectrum,template,z_restframe=3.5,waverange=[8560,8610],plotresult=plotname)

    """
    if verbose: print(' - Loading spectrum and template to cross-correlate ')
    s_wave, s_flux, s_df, s_s2n = felis.load_spectrum(spectrum)

    if waverange is not None:
        goodent = np.where( (s_wave > waverange[0]) & (s_wave < waverange[1]) )
        s_wave, s_flux, s_df, s_s2n = s_wave[goodent], s_flux[goodent], s_df[goodent], s_s2n[goodent]

    if z_restframe is not None:
        s_wave, s_flux, s_df, s_s2n = s_wave / (1+z_restframe), s_flux / (1+z_restframe), s_df, s_s2n
    else:
        z_restframe=1.0

    if spec_median_sub:
        s_flux = s_flux - np.median(s_flux)


    if verbose: print(' - Interpolate template to spectrums wavelength resolution before cross-correlating')
    t_wave_init, t_flux_init, t_df_init, t_s2n_init = felis.load_spectrum(template)
    func       = scipy.interpolate.interp1d(t_wave_init,t_flux_init,kind='linear',fill_value="extrapolate")
    t_flux     = func(s_wave)
    t_flux[t_flux < min_template_level] = 0.0

    if verbose: print(' - Normalizing total template flux to 1')
    if len(t_flux[t_flux != 0]) == 0:
        if verbose: print(' WARNING All interpolated template pixels are 0.0')
    else:
        t_flux     = t_flux / np.sum(t_flux)

    Npix = len(s_flux)

    template_triplelength = np.zeros(3*Npix)
    template_triplelength[0:Npix] = t_flux

    ccresults = np.zeros([Npix,5])

    for ii in np.arange(Npix):
        rollsize       = int(ii+np.floor(Npix/2.))
        template_shift = np.roll(template_triplelength, rollsize)[Npix:-Npix]

        try:
            flux_scale, flux_scale_variance, S2N, chi2_min, NgoodentChi2 = \
                felis.minimize_chi2(s_flux,s_df,template_shift,verbose=False)

            ccresults[ii,:] = flux_scale, flux_scale_variance, S2N, chi2_min, NgoodentChi2
        except:
            print(' ERROR: Problems in minimizing Chi**2 with felis.minimize_chi2() while cross-correlating. '
                  'Stopping for further investigation')

            pdb.set_trace()

        show_illustrative_plot = False
        if show_illustrative_plot:
            xx = np.zeros(24)
            xx[:3] = 1
            plt.plot(xx,'o',markersize=8,color='k')
            plt.plot(np.roll(xx,8)+0.1,'o',markersize=6,color='r')
            plt.plot(np.roll(xx,16)+0.2,'o',markersize=4,color='b')
            plt.plot(np.roll(xx,8)[8:-8]+0.3,'o',markersize=2,color='b')
            plt.plot([7.5,7.5],[0,1.3],'-',color='k')
            plt.plot([15.5,15.5],[0,1.3],'-',color='k')
            plt.show()

    max_S2N     = np.max(ccresults[:,2])
    max_S2N_ent = np.where(ccresults[:,2] == max_S2N)[0][0]
    max_wave    = s_wave[max_S2N_ent]

    # print(' len(temp): '+str(len(t_wave)))
    # print(' len(spec): '+str(len(s_wave)))
    # print(' len(ccresults): '+str(len(ccresults)))
    # print(' CCmac ent: '+str(np.where(ccresults == CCmax)[0]))

    if verbose: print(' - Cross-correlation S2N max happened at index '+str(max_S2N_ent)+
                      ', i.e., lambda = '+str(max_wave)+'A')

    dlam     = np.median(np.diff(s_wave))
    t_wave   = np.arange(np.min(t_wave_init),np.max(t_wave_init),dlam)
    t_cenent = int(np.ceil(len(t_wave)/2))
    max_z    = (max_wave*(z_restframe+1.0) / t_wave[t_cenent]) - 1.0

    template_shift_S2Nmax = np.roll(template_triplelength, int(max_S2N_ent+np.floor(Npix/2.)))[Npix:-Npix]
    flux_scale_S2Nmax     = ccresults[max_S2N_ent,0]

    if plotresult is not None:
        plotname = plotresult
        if verbose: print(' - Setting up and generating plot')
        fig = plt.figure(figsize=(6, 5))
        fig.subplots_adjust(wspace=0.1, hspace=0.5,left=0.1, right=0.99, bottom=0.11, top=0.92)
        Fsize    = 9
        lthick   = 2
        marksize = 4
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif',size=Fsize)
        plt.rc('xtick', labelsize=Fsize)
        plt.rc('ytick', labelsize=Fsize)
        plt.clf()
        plt.ioff()

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        plt.subplot(3,1,1)
        plt.plot(s_wave,template_shift_S2Nmax,'g.',lw=lthick+1, markersize=marksize,alpha=1.0,
                 label='Temp. at z(S/N$_\\textrm{max}$) = '+str(max_z))
        plt.plot(s_wave,flux_scale_S2Nmax * template_shift_S2Nmax,'g-',lw=lthick+2, markersize=marksize,alpha=1.0,
                 label='Temp. flux$_\\textrm{tot}$ scale $\\alpha$ = '+str("%.4f" % flux_scale_S2Nmax)+'')

        plt.plot(s_wave,s_flux,'k-',lw=lthick, markersize=marksize,alpha=1.0,label='Spec.')
        plt.fill_between(s_wave, s_flux-s_df, s_flux+s_df,color='black',alpha=0.2,label='Spec. err')

        plt.plot([max_wave,max_wave],[np.min(s_flux-s_df),np.max(s_flux+s_df)],'--r',lw=lthick,
                 markersize=marksize,alpha=1.0,label='S/N$_\\textrm{max}$ = '+str("%.4f" % max_S2N)+'')
        plt.xlabel(' Wavelength [A]')
        plt.ylabel(' Flux ')
        #plt.ylim([-2,6])
        leg = plt.legend(fancybox=True, loc='upper center',prop={'size':Fsize/1.3},ncol=3,numpoints=1,
                         bbox_to_anchor=(0.5, 1.43))  # add the legend
        leg.get_frame().set_alpha(0.7)

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        plt.subplot(3,1,2)
        plt.plot(s_wave,ccresults[:,2],'-r',lw=lthick, markersize=marksize,alpha=1.0)
        plt.plot([max_wave,max_wave],[np.min(ccresults[:,2]),np.max(ccresults[:,2])],'--r',lw=lthick,
                 markersize=marksize,alpha=1.0)
        plt.xlabel(' Wavelength [A]')
        plt.ylabel(' Cross-Correlation S/N')
        # leg = plt.legend(fancybox=True, loc='upper right',prop={'size':Fsize},ncol=1,numpoints=1)
        # leg.get_frame().set_alpha(0.7)

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        plt.subplot(3,1,3)
        plt.plot(s_wave,ccresults[:,3],'-r',lw=lthick, markersize=marksize,alpha=1.0)
        plt.plot([max_wave,max_wave],[np.min(ccresults[:,3]),np.max(ccresults[:,3])],'--r',lw=lthick,
                 markersize=marksize,alpha=1.0)
        plt.xlabel(' Wavelength [A]')
        plt.ylabel(' $\chi^2_\\textrm{min}$')
        # plt.ylim([0,1])
        # leg = plt.legend(fancybox=True, loc='upper right',prop={'size':Fsize},ncol=1,numpoints=1)
        # leg.get_frame().set_alpha(0.7)

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        if verbose: print('   Saving plot to '+plotname)
        plt.savefig(plotname)
        plt.clf()
        plt.close('all')

    return  s_wave, ccresults, max_S2N, max_z
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def cross_correlate_wscipy(spectrum,template,z_restframe=None,ccS2N=False,waverange=None,plotresult=False,verbose=True):
    """
    Cross-correlate a spectrum with a template.
    To handle wavelength resolution properly, move to rest-frame (otherwise line-spacings etc.
    will likely be off for cross-correlation).

    For the rest-frame case the template is interpolated to ensure such an agreement.

    NB! This cross-correlation function is superseded by 'match_template' which does a manual template matching
        and cross-correlation using minimize_chi2() and cross_correlate_template()

    --- INPUT ---
    spectrum       fits spectrum to find a (cross-correlation) match to template for
    template       fits template to correlate with
    z_restframe    To perform cross-correlation in rest-frame (shifting the spectrum) provide a redshift.
    ccS2N          Use cross-correlation S/N to determine maximum value and redshift of peak instead of flux?
                   Both CC flux and variance will be returned for post analysis.
    waverange      To prevent running the cross-correlation on the full spectrum,
                   provide the wavelength range to restrict the correlation to.
    plotresult     Provide an output file name to illustrate the location for the best
                   correlation between the spectrum and template
    verbose        Toggle verbosity

    --- EXAMPLE OF USE ---
    import felis

    specdir  = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/felis_testing/'
    spectrum = specdir+'uves_felis_mock_MUSEspectrum_noisesigma1p0.fits'
    template = specdir+'uves_felis_template_CIIIdoublet_sig_0p5_fluxCIII1_2p0_fluxratio_0p5.fits'
    #template = specdir+'uves_felis_template_singleline_Singleline_sig_0p5_flux_10p0.fits'

    plotname = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/felis_testing/test_CCResultPlot.pdf'

    wave, cc_flux, cc_variance, cc_norm, z_CCmax  = felis.cross_correlate(spectrum,template,waverange=[8500,8650],plotresult=plotname,ccS2N=True)

    check normalization with: scipy.integrate.trapz(ccresults, x=wave)

    """
    if verbose: print(' - Loading spectrum and template to cross-correlate ')
    s_wave, s_flux, s_df, s_s2n = felis.load_spectrum(spectrum)

    if waverange is not None:
        goodent = np.where( (s_wave > waverange[0]) & (s_wave < waverange[1]) )
        s_wave, s_flux, s_df, s_s2n = s_wave[goodent], s_flux[goodent], s_df[goodent], s_s2n[goodent]

    if z_restframe is not None:
        s_wave, s_flux, s_df, s_s2n = s_wave / (1+z_restframe), s_flux / (1+z_restframe), s_df, s_s2n
    else:
        z_restframe=1.0

    if verbose: print(' - Interpolate template to spectrums wavelength resolution before cross-correlating')
    t_wave_init, t_flux_init, t_df_init, t_s2n_init = felis.load_spectrum(template)
    func       = scipy.interpolate.interp1d(t_wave_init,t_flux_init,kind='linear',fill_value=0)
    dlam       = np.median(np.diff(s_wave))
    t_wave     = np.arange(np.min(t_wave_init),np.max(t_wave_init),dlam)
    t_flux     = func(t_wave)

    print(template+'\n ------> Max(template) = '+str(np.max(t_flux)))
    try:
        ccresults_flux      = scipy.signal.correlate(s_flux,  t_flux,    mode='same', method='auto') # cc flux
        ccresults_variance  = scipy.signal.correlate(s_df**2, t_flux**2, mode='same', method='auto') # cc noise
        # scipy.signal.correlate = scipy.signal.convolve are identical when template is symmetric.
        # Otherwise convolution convolves (smoothes) the signal by the filter (i.e., tries to answer:
        # "what is the output of this filter if the input is x(t)"), whereas the cross-correlations
        # tries to answer: "Given a spectrum, is the signal represented by the template, somehow present in
        # the spectrum?"
        # https://dsp.stackexchange.com/questions/27451/the-difference-between-convolution-and-cross-correlation-from-a-signal-analysis
        # https://www.youtube.com/watch?v=C3EEy8adxvc
    except:
        print(' WARNING Working on '+spectrum)
        print('         but ran into problems with scipy.signal.correlate; stopping for investigation')
        pdb.set_trace()

    s_flux_norm = (s_flux - np.mean(s_flux)) / (np.std(s_flux) * len(s_flux))
    t_flux_norm = (t_flux - np.mean(t_flux)) /  np.std(t_flux)
    ccresults_norm = scipy.signal.correlate(s_flux_norm, t_flux_norm, mode='same', method='auto')
    #ccresults_norm = ccresults_norm / scipy.integrate.trapz(ccresults_norm, x=s_wave)


    # Correlation method described at
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.correlate.html#scipy.signal.correlate
    CCmax_flux     = np.max(ccresults_flux)
    CCmax_S2N      = np.max(ccresults_flux/np.sqrt(ccresults_variance))

    if ccS2N:
        if verbose: print(' - Choosing cross-correlation max and signal: Using CC_S/N expressions')
        cctype         = 'S/N'
        CCmax          = CCmax_S2N
        ccresults      = ccresults_flux / np.sqrt(ccresults_variance)
    else:
        if verbose: print(' - Choosing cross-correlation max and signal: Using CC_flux expressions')
        cctype         = 'flux'
        CCmax          = CCmax_flux
        ccresults      = ccresults_flux

    try:
        if len(CCmax) > 1:
            print(' WARNING: Cross-correlation has a CC_S/N maximum at '+str(len(CCmax))+' wavelengths.')
            print('          Will use first (index '+str(CCmax[0])+') for z_CCmax estimate ')
            CCmax = CCmax[0]
    except:
        pass

    if len(ccresults[np.isfinite(ccresults)]) == 0:
        CCmax_ent = 0
    else:
        CCmax_ent = np.where(ccresults == CCmax)[0]
        CCmax_ent = CCmax_ent[0]
        # print(' len(temp): '+str(len(t_wave)))
        # print(' len(spec): '+str(len(s_wave)))
        # print(' len(ccresults): '+str(len(ccresults)))
        # print(' CCmac ent: '+str(np.where(ccresults == CCmax)[0]))

    if verbose: print(' - Cross-correlation max happened at index '+str(CCmax_ent)+', i.e., lambda = '+str(s_wave[CCmax_ent])+'A')

    t_cenent = int(np.ceil(len(t_wave)/2))
    z_CCmax  = (s_wave[CCmax_ent]*(z_restframe+1.0) / t_wave[t_cenent]) - 1.0

    if plotresult is not None:
        plotname = plotresult
        if verbose: print(' - Setting up and generating plot')
        fig = plt.figure(figsize=(7, 6))
        fig.subplots_adjust(wspace=0.1, hspace=0.5,left=0.1, right=0.97, bottom=0.11, top=0.92)
        Fsize    = 12
        lthick   = 2
        marksize = 4
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif',size=Fsize)
        plt.rc('xtick', labelsize=Fsize)
        plt.rc('ytick', labelsize=Fsize)
        plt.clf()
        plt.ioff()

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        plt.subplot(3,1,1)
        waveshift = t_wave - t_wave[t_cenent] + s_wave[CCmax_ent]
        plt.plot(waveshift,t_flux,'g-',lw=lthick+2, markersize=marksize,alpha=1.0,
                 label='Template at z(CCmax) = '+str(z_CCmax))
        plt.plot(s_wave,s_flux,'k-',lw=lthick, markersize=marksize,alpha=1.0,label='Spectrum')
        plt.plot([s_wave[CCmax_ent],s_wave[CCmax_ent]],[np.min(s_flux),np.max(s_flux)],'--r',lw=lthick,
                 markersize=marksize,alpha=1.0,label='Max(CC Value)')
        plt.xlabel(' Wavelength [A]')
        plt.ylabel(' Flux ')
        leg = plt.legend(fancybox=True, loc='upper center',prop={'size':Fsize/1.3},ncol=3,numpoints=1,
                         bbox_to_anchor=(0.5, 1.3))  # add the legend
        leg.get_frame().set_alpha(0.7)

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        plt.subplot(3,1,2)
        plt.plot(s_wave,ccresults,'-r',lw=lthick, markersize=marksize,alpha=1.0)
        plt.plot([s_wave[CCmax_ent],s_wave[CCmax_ent]],[np.min(ccresults),np.max(ccresults)],'--r',lw=lthick,
                 markersize=marksize,alpha=1.0,label='Max(CC Value)')
        plt.xlabel(' Wavelength [A]')
        plt.ylabel(' Linear CC ('+cctype+')')
        # leg = plt.legend(fancybox=True, loc='upper right',prop={'size':Fsize},ncol=1,numpoints=1)
        # leg.get_frame().set_alpha(0.7)

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        plt.subplot(3,1,3)
        plt.plot(s_wave,ccresults_norm,'-r',lw=lthick, markersize=marksize,alpha=1.0)
        plt.plot([s_wave[CCmax_ent],s_wave[CCmax_ent]],[np.min(ccresults_norm),np.max(ccresults_norm)],'--r',lw=lthick,
                 markersize=marksize,alpha=1.0,label='Max(CC Value)')
        plt.xlabel(' Wavelength [A]')
        plt.ylabel(' Linear CC (flux) Normalized')
        plt.ylim()
        # leg = plt.legend(fancybox=True, loc='upper right',prop={'size':Fsize},ncol=1,numpoints=1)
        # leg.get_frame().set_alpha(0.7)

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        if verbose: print('   Saving plot to '+plotname)
        plt.savefig(plotname)
        plt.clf()
        plt.close('all')

    return  s_wave, ccresults_flux, ccresults_variance, ccresults_norm, z_CCmax
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def load_spectrum(specfile,verbose=True):
    """
    Loading spectrum generated with felis.save_spectrum() or TDOSE.

    --- INPUT ---

    --- EXAMPLE OF USE ---
    import felis
    specfile  = './spectrum_output.fits'
    w, f, df, s2n = felis.load_spectrum(specfile)

    """
    if verbose: print(' - Loading SPEC1D extension (spectrum) of \n   '+specfile)
    dat = fits.open(specfile)['SPEC1D'].data

    wave    = dat['wave']
    flux    = dat['flux']
    fluxerr = dat['fluxerror']
    s2n     = dat['s2n']

    return wave, flux, fluxerr, s2n
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
def save_dictionary(dictionary, output='./saveddictionary_RENAME_.pkl'):
    with open(output, 'wb') as f:
        pickle.dump(dictionary, f, pickle.HIGHEST_PROTOCOL)
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def load_picklefile(picklefile):
    with open(picklefile, 'rb') as f:
        return pickle.load(f)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
