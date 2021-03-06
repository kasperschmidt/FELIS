# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
import felis
from astropy.io import fits
import scipy.signal
import pdb
import matplotlib.pyplot as plt
import datetime
import pickle
import sys
import os
import numpy as np
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def match_templates2specs(templates,spectra,speczs,picklename,wavewindow=[50.0],wavecen_restframe=[1908.0],
                          vshift=None,min_template_level=1e-4,plotdir=None,plot_allCCresults=False,
                          subtract_spec_median=True,overwrite=False,verbose=True):
    """
    Wrapper around felis cross-correlation template matching, to match a list of spectra with a list of templtes.

    --- INPUT ---
    spectra               fits spectra to find a (cross-correlation) match to template for
    templates             fits templates to correlate with
    speczs                Spectroscopic redshifts to perform cross-correlation in rest-frame (shifting the spectrum).
    subtract_spec_median  Subtract median value of spectrum (approximating the continuum level)
    picklename            Name of pickle file to store final cross-correlation results in
    wavewindow            Window (wavecen_restframe * (1+speczs) +/- wavewindow) to perform template matching over.
    wavecen_restframe     Central rest-frame  wavelength of the region to match
    vshift                If a velcotiy shift is known, provide it here and it will be stored in output (not used)
    min_template_level    The template is interpolated to the wavelength grid of the spectrum and extrapolated
                          beyond it's edges if nescessary. In this extrapolation (assuming the template goes to ~0
                          at the edges), very small values (e.g., <1e-20) can be returned. To set these to 0.0
                          provide a level below which all values in the interpolated template are treated as 0s.
    plotdir               Directory to store plots to
    plot_allCCresults     To plot all the cross-correlation plots, set this to True
    overwrite             Overwrite existing pickle file if it already exists?
    verbose               Toggle verbosity

    --- EXAMPLE OF USE ---
    import felis
    import glob

    specdir  = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/felis_testing/'
    #specs    = glob.glob(specdir+'uves_felis_mock_MUSEspectrum_noisesigma*3p0.fits')
    specs    = glob.glob(specdir+'uves_felis_mock_MUSEspectrum_noisesigma*.fits')
    speczs   = [3.5]*len(specs)

    tempdir  = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/felis_testing/'
    #temps    = glob.glob(specdir+'uves_felis_template_CIIIdoublet_sig_0p25_fluxCIII1_4p0_flux*.fits')
    temps    = glob.glob(specdir+'uves_felis_template_CIIIdoublet_*fits')
    temps    = glob.glob(specdir+'uves_felis_template_CIVdoublet_*fits')

    plotdir  = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/felis_testing/plots_CCresults180615/'
    pickle   = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/felis_testing/CCresults180615_RENAME_.pkl'

    specs = ['/Volumes/DATABCKUP1/TDOSEextractions/171201_TDOSEextraction/Modelimg/tdose_spectra/tdose_spectrum_candels-cdfs-15_modelimg_0115003085-0115003085.fits']
    speczs = [3.2585198879241943]

    ccdic    = felis.match_templates2specs(temps,specs,speczs,pickle,wavewindow=[60]*len(specs),plotdir=plotdir,wavecen_restframe=[1549.0]*len(specs))

    --- OUTPUT ---

    This wrapper will collect all the cross-correlation results in a main dictionary.
    The dictionary will be returned directly but also saved to disk as a pickled filed.
    This file can be loaded with: felis.load_picklefile(picklefilename)

    The returned dictionary has the following format:

    dictionary.keys()  = the list of spectra that have been crossmatched (input = 'spectra')

    Each entry in the dictionary (dictionary[key]) contains the following entries:
        'wavelengths'               :   The wavelength vector used for the cross-correlation matching of
                                        each of the N templates
        'templatevec'               :   list of each of the N templates matched to spectrum
        'zspec'                     :   Spectroscopic redshift for spectrum
        'zCCmaxvec'                 :   the redshift corresponding to max S/N for each of the N templates matched
        'ccresultsarray_flux'       :   fluc vectors for each of the N templates matched
        'ccresultsarray_variance'   :   Variance vecotrs for each of the N templates matched
        'ccresultsarr_S2N'          :   S/N vector for each of the N templates matched
        'ccresultsarr_chi2'         :   chi squared values for the cross-correlation of each of the N templates matched
        'ccresultsarr_Ngoodent'     :   The number of good pixels used in the cross correlation for each of
                                        the N templates matched
        'S2NCCmaxvec'               :   vector with max(S/N) values for N templates matched
        'continuumlevel'            :   The 'continuum level' of the spectrum removed in the cross-correlation.
                                        Currently the value is simply the median of the spectrum.
        'vshift'                    :   If a velocity shift was provided for the template match this is stored here

    The picklefile can be used to assemble sub-sample results (e.g., S/N cuts) based on the
    template cross-correlations with
        felis.selection_from_picklefile()
    And individual entries can be plotted using
        felis.plot_picklefilecontent()

    """
    ccresultdic = {}

    if verbose: print(' - Starting cross-correlation of the '+str(len(spectra))+' spectra and '+
                      str(len(templates))+' templates')
    startstring = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S.%f")
    if verbose: print('   '+startstring+'\n')

    if len(spectra) == 0:
        sys.exit(' No spectra provided')

    if len(templates) == 0:
        sys.exit(' No templates provided')

    for ss, spec in enumerate(spectra):
        # Nwave = pyfits.open(spec)[1].header['NAXIS2']
        spec_namebase = spec.split('/')[-1].split('.fit')[0]
        for tt, temp in enumerate(templates):
            temp_namebase = temp.split('/')[-1].split('.fit')[0]

            wavecenter  = wavecen_restframe[ss]  * (1.0 + speczs[ss])
            waverange   = [wavecenter-wavewindow[ss],wavecenter+wavewindow[ss]]

            wave, ccresults, max_S2N, max_z, continuumlevel = \
                felis.cross_correlate_template(spec,temp,z_restframe=speczs[ss],spec_median_sub=subtract_spec_median,
                                               waverange=waverange,min_template_level=min_template_level,verbose=verbose)

            if tt == 0:
                ccresultsarr_flux      = np.array(ccresults[:,0])
                ccresultsarr_variance  = ccresults[:,1]
                ccresultsarr_S2N       = ccresults[:,2]
                ccresultsarr_chi2      = ccresults[:,3]
                ccresultsarr_Ngoodent  = ccresults[:,4]
                templatevec            = np.array([temp])
                S2NCCmaxvec            = np.array([max_S2N])
                zCCmaxvec              = np.array([max_z])
                #print('-------------------------->'+str(np.max(ccresultsarr_S2N))+'  '+str(max_S2N))
            else:
                ccresultsarr_flux      = np.vstack((ccresultsarr_flux,ccresults[:,0]))
                ccresultsarr_variance  = np.vstack((ccresultsarr_variance,ccresults[:,1]))
                ccresultsarr_S2N       = np.vstack((ccresultsarr_S2N,ccresults[:,2]))
                ccresultsarr_chi2      = np.vstack((ccresultsarr_chi2,ccresults[:,3]))
                ccresultsarr_Ngoodent  = np.vstack((ccresultsarr_Ngoodent,ccresults[:,4]))
                templatevec            = np.append(templatevec,temp)
                S2NCCmaxvec            = np.append(S2NCCmaxvec,max_S2N)
                zCCmaxvec              = np.append(zCCmaxvec,max_z)
                #print('-------------------------->'+str(np.max(ccresultsarr_S2N[tt,:]))+'  '+str(S2NCCmaxvec[tt]))

        ccresultdic[spec]  = {'wavelengths':wave, 'templatevec':templatevec, 'zspec':speczs[ss],
                              'zCCmaxvec':zCCmaxvec, 'ccresultsarray_flux':ccresultsarr_flux,
                              'ccresultsarray_variance':ccresultsarr_variance,'S2NCCmaxvec':S2NCCmaxvec,
                              'ccresultsarr_S2N':ccresultsarr_S2N,'ccresultsarr_chi2':ccresultsarr_chi2,
                              'ccresultsarr_Ngoodent':ccresultsarr_Ngoodent}

        ccresultdic[spec]['continuumlevel'] = continuumlevel
        if vshift is not None:
            ccresultdic[spec]['vshift'] = vshift[ss]


    if verbose: print('\n - Finished cross-correlation of the '+str(len(spectra))+' spectra and '+
                      str(len(templates))+' templates')
    if verbose: print('   '+datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S.%f")+'  (started at '+startstring+')')

    if verbose: print(' - Saving dictionary to '+picklename)
    if os.path.isfile(picklename) & (overwrite == False):
        print('\n   FELIS WARNING: The output file '+picklename+'exists but overwrite==False so returning "None" ')
        return None
    else:
        felis.save_dictionary(ccresultdic,picklename)

        if verbose: print(' - Will plot all cross-correlation results in saved dictionary as plot_allCCresults=True')
        if plot_allCCresults:
            for ss, spec in enumerate(spectra):
                spec_namebase = spec.split('/')[-1].split('.fit')[0]
                for tt, temp in enumerate(templates):
                    temp_namebase = temp.split('/')[-1].split('.fit')[0]
                    plotname      = plotdir+spec_namebase+'_CCwith_'+temp_namebase+'.pdf'

                    felis.plot_picklefilecontent([spec],picklename,plottemplates=[tt],showspecerr=False,
                                                 plotnames=[plotname],plotdir=plotdir,verbose=verbose)

        if verbose: print(' - Returning dictioniary with results ')
        loaddic = felis.load_picklefile(picklename)
        return  loaddic
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def selection_from_picklefile(picklefile,S2Nmaxrange=[3,1e4],zspecrange=[0,10],voffsetrange=[-1e4,1e4],
                              zspecISzLya=False,verbose=True):
    """
    Function, returning the list of spectra (keys) from the FELIS tempalte match pickled dictionary satisfying
    a set of criteria on SN, redshift, velocity offset, etc.

    --- INPUT ---
    picklefile          The path and name to FELIS output pickelfile to select subset of matches from.
    S2Nmaxrange         The range of S/N template accepted for template matches returned.
    zspecrange          The range of redshifts accepted for template matches returned.
    voffsetrange        The range of velocity offsets in km/s (wrt. to the spectroscopic redshift in dictionary)
                        accepted for template matches returned.
    zspecISzLya         If picklefile contains zLya keyword insead of zspec (was generated before 180912), set this
                        keyword to True to enable proper handling of the dictionary keywords.
    verbose             Toggle the verbosity.

    --- EXAMPLE OF USE ---
    import felis

    picklepath = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/'
    picklefile = picklepath+'MUSEWideLAEs_CCresults180502_CIV_all_575specX12templates.pkl'

    goodkeys   = felis.selection_from_picklefile(picklefile,S2Nmaxrange=[3,5])

    """
    if verbose: print(' - Loading the picklefile \n   '+picklefile)
    loaddic  = felis.load_picklefile(picklefile)
    goodkeys = []
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Ensure compatibility with FELIS output dictionaries from before 180912
    zkey = 'zspec'
    if zspecISzLya: zkey = 'zLya'
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Looking for keys in the pickled dictionary satisfying the cuts:')
    if verbose: print('       '+zkey+'          :  ['+str(zspecrange[0])+','+str(zspecrange[1])+']')
    if verbose: print('       max(S/N)       :  ['+str(S2Nmaxrange[0])+','+str(S2Nmaxrange[1])+']')
    if verbose: print('       voffset[km/s]  :  ['+str(voffsetrange[0])+','+str(voffsetrange[1])+']')
    for key in loaddic.keys():
        keydic = loaddic[key]

        template, vshift_intr, vshift_match, flux, fluxerr, S2Nmax, Ngoodent, chi2, zspec, zS2Nmax =  \
            felis.getresult4maxS2N(loaddic,key,zspecISzLya=zspecISzLya)

        if ((keydic[zkey] > zspecrange[0]) & (keydic[zkey] < zspecrange[1])) & \
                ((S2Nmax > S2Nmaxrange[0]) & (S2Nmax < S2Nmaxrange[1])) & \
                ((vshift_match > voffsetrange[0]) & (vshift_match < voffsetrange[1])):
            goodkeys.append(key)
    if verbose: print(' - Found '+str(len(goodkeys))+' keys in the pickled dictionary satisfying the cuts:')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Returning those')
    return goodkeys
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_picklefilecontent(specs2plot,picklefile,plotnames=None,plotdir=None,z_restframe=None,
                           plottemplates=None,showspecerr=True,zspecISzLya=False,verbose=True):
    """
    Function to plot individual results from a picklefile generated with felis.match_templates2specs()

    --- INPUT ---
    specs2plot          The spectra from the pickle file (the pickle dictinary keys) to plot. These can be
                        provide in a list by hand or be selected with felis.selection_from_picklefile()
    picklefile          The path and name to pickelfile to plot content of
    plotnames           The names of the plot(s) to generate. If 'None' the string '_templatematch.pdf' will
                        be appended that pickle dictionary key (spectrum name).
    plotdir             To chose a different directory for saving the plots (other than the directory
                        in which the spectrum provided is stored) provide this here.
    z_restframe         The redshift used to move spectrum to rest-frame
    plottemplates       By default the template match with highest S/N is plotted. To plot another template
                        provide the entry in the dictionary of this template. To see the templates, look at
                        the 'templatevec' entry:
                            picload = felis.load_picklefile(picklefile)
                            picload[specname]['templatevec']
    showspecerr         Show the error on the data spectrum? Can make the automatically set y-axis range less ideal.
    zspecISzLya         If picklefile contains zLya keyword insead of zspec (was generated before 180912), set this
                        keyword to True to enable proper handling of the dictionary keywords.
    verbose             Toggle the verbosity.

    --- EXAMPLE OF USE ---
    import felis

    picklepath = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/'
    picklefile = picklepath+'MUSEWideLAEs_CCresults180502_CIV_all_575specX12templates.pkl'

    picklefile = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/felis_testing/CCresults180615_RENAME_.pkl'

    plotdir    = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/MUSEwideLAE_FELISplots/'
    specs2plot = ['/Volumes/DATABCKUP1/TDOSEextractions/171201_TDOSEextraction/Modelimg/tdose_spectra/tdose_spectrum_candels-cdfs-15_modelimg_0115003085-0115003085.fits']

    felis.plot_picklefilecontent(specs2plot,picklefile,plotnames=None,plotdir=plotdir,verbose=True)


    """
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Ensure compatibility with FELIS output dictionaries from before 180912
    zkey = 'zspec'
    if zspecISzLya: zkey = 'zLya'
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - loading the dictionary from the pickled file:\n   '+picklefile)
    loaddic  = felis.load_picklefile(picklefile)

    for ss, spec in enumerate(specs2plot):
        if plotnames is None:
            plotname = spec.replace('.fits','_templatematch.pdf')
        else:
            plotname = plotnames[ss]

        if plotdir is not None:
            specdir  = '/'.join(plotname.split('/')[:-1])
            plotname = plotname.replace(specdir,plotdir)

        spec_wave, spec_flux, spec_df, spec_s2n = felis.load_spectrum(spec,verbose=verbose)
        spec_dic  = loaddic[spec]

        if plottemplates is None:
            # Getting the entry of the template with maximum S/N in the CC results
            max_S2N = np.max(spec_dic['S2NCCmaxvec'])
            besttemplate_ent = np.where(spec_dic['ccresultsarr_S2N'] == max_S2N)[0][0] # Template of max S/N
        else:
            max_S2N          = spec_dic['S2NCCmaxvec'][plottemplates[ss]]
            besttemplate_ent = plottemplates[ss]

        # moving spectrum to restframe
        if z_restframe is None:
            z_spec=spec_dic[zkey]

            spec_wave, spec_flux, spec_df, spec_s2n = \
                spec_wave / (1+z_spec), spec_flux * (1.0+z_spec), spec_df * (1.0+z_spec), spec_s2n
        else:
            z_spec = 0.0

        # subtract continuum level from spectrum
        spec_flux = spec_flux - spec_dic['continuumlevel'] * (1.0+z_spec)

        # Limit spectrum to range of wavelengths cross-correlated with tempalte
        goodent = np.where( (spec_wave >= spec_dic['wavelengths'][0]) & (spec_wave <= spec_dic['wavelengths'][-1]) )[0]
        spec_wave, spec_flux, spec_df, spec_s2n = \
            spec_wave[goodent], spec_flux[goodent], spec_df[goodent], spec_s2n[goodent]

        template         = spec_dic['templatevec'][besttemplate_ent]
        max_z            = spec_dic['zCCmaxvec'][besttemplate_ent]

        t_wave_init, t_flux_init, t_df_init, t_s2n_init = felis.load_spectrum(template,verbose=verbose)
        func       = scipy.interpolate.interp1d(t_wave_init,t_flux_init,kind='linear',fill_value="extrapolate")
        t_flux     = func(spec_wave)

        # Getting the entry in the CC flux scalings vector for the given template where S/N is max
        max_S2N_ent      = np.where(spec_dic['ccresultsarr_S2N'][besttemplate_ent,:] == max_S2N)[0][0]

        Npix = len(spec_flux)
        template_triplelength            = np.zeros(3*Npix)
        template_triplelength[0:Npix]    = t_flux
        template_shift_S2Nmax            = np.roll(template_triplelength, int(max_S2N_ent+np.floor(Npix/2.)))[Npix:-Npix]
        template_shift_S2Nmax_normalized = template_shift_S2Nmax/np.trapz(template_shift_S2Nmax,spec_wave)

        flux_scale_S2Nmax     = spec_dic['ccresultsarray_flux'][besttemplate_ent,max_S2N_ent]
        max_wave              = spec_dic['wavelengths'][max_S2N_ent]

        if verbose: print(' - Setting up and generating plot:\n   '+plotname)
        fig = plt.figure(figsize=(6, 7))
        fig.subplots_adjust(wspace=0.1, hspace=0.5,left=0.1, right=0.99, bottom=0.07, top=0.91)
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
        ax = plt.subplot(4,1,1)

        plt.plot(spec_wave,template_shift_S2Nmax_normalized,'g.',lw=lthick+1, markersize=marksize,alpha=1.0,
                 label='Temp. at z(S/N$_\\textrm{max}$) = '+str(max_z))
        plt.plot(spec_wave,flux_scale_S2Nmax * template_shift_S2Nmax_normalized,
                 'g-',lw=lthick+2, markersize=marksize,alpha=1.0,
                 label='Temp. flux$_\\textrm{tot}$ scale $\\alpha$ = '+str("%.4f" % flux_scale_S2Nmax)+'')

        Ftot_trapz = np.trapz(flux_scale_S2Nmax * template_shift_S2Nmax_normalized,spec_wave)

        if spec_dic['continuumlevel'] != 0.0:
            plt.plot(spec_wave,spec_flux+spec_dic['continuumlevel']*(1+z_spec),'k-',lw=lthick,
                     markersize=marksize,alpha=0.5, label='Spectrum (w/ cont.)')
            speclabel = 'Spectrum (w/o cont.)'
        else:
            speclabel = 'Spectrum (w/ cont.)'
        plt.plot(spec_wave,spec_flux,'k-',lw=lthick, markersize=marksize,alpha=1.0,label=speclabel)

        #what about...
        #  plt.step(specdat['wave'],specdat['flux'],'r',where='mid')
        #  plt.fill_between(specdat['wave'],specdat['flux']-specdat['fluxerror'],specdat['flux']+specdat['fluxerror'],step='mid',color='red',alpha=0.5)



        if showspecerr:
            plt.fill_between(spec_wave, spec_flux-spec_df, spec_flux+spec_df,color='black',alpha=0.2,label='Spectrum err')

            SNlineYrange = [np.min(spec_flux-spec_df),np.max(spec_flux+spec_df)]
        else:
            SNlineYrange = [np.min(spec_flux),np.max(spec_flux)]

        plt.plot([max_wave,max_wave],SNlineYrange,'--r',lw=lthick,
                 markersize=marksize,alpha=1.0,label='S/N$_\\textrm{max}$ = '+str("%.4f" % max_S2N)+'')

        plt.xlabel(' Wavelength [A]')
        plt.ylabel(' Flux ')
        #plt.ylim([-2,6])
        leg = plt.legend(fancybox=True, loc='upper center',prop={'size':Fsize/1.3},ncol=3,numpoints=1,
                         bbox_to_anchor=(0.45, 1.43))  # add the legend
        leg.get_frame().set_alpha(0.7)
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        plt.text(0.45, 1.45, 'Template No. '+str(int(besttemplate_ent+1))+'/'+str(len(spec_dic['templatevec']))+
                 ': '+template.split('/')[-1].replace('_','\_')+
                 ' (Ftot\_trapz = '+str("%.2f" % Ftot_trapz)+')',
                 fontsize=Fsize/1.3, horizontalalignment='center', transform=ax.transAxes)

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        plt.subplot(4,1,2)

        plt.plot(spec_wave,spec_s2n,'k-',lw=lthick, markersize=marksize,alpha=0.5)
        plt.plot(spec_wave,spec_flux/spec_df,'k-',lw=lthick, markersize=marksize,alpha=1.0)

        plt.plot([max_wave,max_wave],[np.min(spec_flux/spec_df),np.max(spec_flux/spec_df)],'--r',lw=lthick,
                 markersize=marksize,alpha=1.0)
        plt.xlabel(' Wavelength [A]')
        plt.ylabel(' Spectrum S/N ')

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        plt.subplot(4,1,3)
        plt.plot(spec_dic['wavelengths'],spec_dic['ccresultsarr_S2N'][besttemplate_ent,:],'-r',lw=lthick, markersize=marksize,alpha=1.0)
        plt.plot([max_wave,max_wave],[np.min(spec_dic['ccresultsarr_S2N'][besttemplate_ent,:]),
                                      np.max(spec_dic['ccresultsarr_S2N'])],'--r',lw=lthick,
                 markersize=marksize,alpha=1.0)
        plt.xlabel(' Wavelength [A]')
        plt.ylabel(' Cross-Correlation S/N')

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        plt.subplot(4,1,4)
        plt.plot(spec_dic['wavelengths'],spec_dic['ccresultsarr_chi2'][besttemplate_ent,:],'-r',lw=lthick, markersize=marksize,alpha=1.0)
        plt.plot([max_wave,max_wave],[np.min(spec_dic['ccresultsarr_chi2'][besttemplate_ent,:]),np.max(spec_dic['ccresultsarr_chi2'][besttemplate_ent,:])],'--r',lw=lthick,
                 markersize=marksize,alpha=1.0)
        plt.xlabel(' Wavelength [A]')
        plt.ylabel(' $\chi^2_\\textrm{min}$')

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        if verbose: print('   Saving plot to '+plotname)
        plt.savefig(plotname)
        plt.clf()
        plt.close('all')

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def minimize_chi2(data,data_sigma,template,temp_wave=None,verbose=True):
    """
    Minimizing chi2 for a 1D template matching and returning the corresponding (flux) scaling, alpha.

    Done using that
        chi**2         = Sum_i (d_i - alpha * t_i)**2 / sigma_i**2
    and that dchi**2/dalpha = 0 implies
        alpha          = ( Sum_i d_i*t_i/sigma_i**2 ) / ( Sum_i t_i**2 / sigma_i**2 )
        sigma_alpha**2 = 1 / Sum_ix( t_i**2 / sigma_i**2 ) = alpha_variance

    Here, d is the data (flux from spectrum), sigma is the uncertainty on the data, t is the template,
    alpha is the flux scaling of the template, sigma_alpha is the uncertainty on alpha, and i runs
    over the pixels in the spectrum.

    --- INPUT ---

    data            Data to match template to
    data_sigma      Uncertainty on data, i.e., sqrt(variance)
    template        Normalized template to search for in data. Should be of the same lenght as data
    temp_wave       Provide template wavelengths for normalization check if dlam is not 1.
    verbose         Toggle verbosity

    """
    normlim = 1e-10
    if temp_wave is not None:
        temp_int = np.trapz(template,temp_wave)
    else:
        temp_int = np.trapz(template)

    if temp_int == 0: # all pixels are 0
        pass
    else:
        sumdiff = np.abs(temp_int-1.0)
        if sumdiff > normlim: #checking if provided template is normalized
            raise ValueError('FELIS WARNING: Template is not normalized: |np.trapz(template,dwave)-1.0| = '+
                             str(sumdiff)+' > '+str(normlim))

    if len(data) == len(template):
        Npix = len(data)
    else:
        sys.exit('The length of the data and template should be the same; it is not.')
    if verbose: print(' - Will calculate and minimize chi**2 for data and template of length '+str(Npix))
    goodent  = np.where((data_sigma > 0) & (np.isfinite(data)) & (template != 0))[0]
    Ngoodent = len(goodent)

    if Ngoodent == 0:
        if verbose: print(' - No entries left where data_sigma > 0 & data is finite & template != 0, so returning 0s')
        chi2_min        = 0.0
        alpha           = 0.0
        alpha_variance  = 0.0
        S2N             = 0.0
    else:
        t2err2         = template[goodent]**2 / data_sigma[goodent]**2
        alpha          = np.sum( data[goodent] * template[goodent] / data_sigma[goodent]**2 ) / np.sum(t2err2)
        alpha_variance = 1.0 / np.sum(t2err2)

        S2N            = alpha / np.sqrt(alpha_variance)
        chi2_min       = np.sum( ( data[goodent]-alpha*template[goodent] )**2 / data_sigma[goodent]**2 )

    return alpha, alpha_variance, S2N, chi2_min, Ngoodent
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def cross_correlate_template(spectrum,template,z_restframe=None,waverange=None,
                             min_template_level=1e-2,spec_median_sub=False,verbose=True):
    """
    Function to cross-correlate a spectrum with a template using felis.minimize_chi2().
    It can be run for multiple spectra with the wrapper match_templates2specs which also enables
    easy plotting of the results.

    --- INPUT ---
    spectrum            fits spectrum to find a (cross-correlation) match to template for (flux in f_lambda)
    template            fits template to correlate with
    z_restframe         To perform cross-correlation in rest-frame (shifting the spectrum) provide a redshift.
    waverange           To prevent running the cross-correlation on the full spectrum, provide the wavelength
                        range (based on wavelength-vector in spectrum) to restrict the correlation to.
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

    s_wave, ccresults, max_S2N, max_z, continuumlevel = felis.cross_correlate_template(spectrum,template,z_restframe=3.5,waverange=[8560,8610])

    """
    if verbose: print(' - Loading spectrum and template to cross-correlate ')
    s_wave, s_flux, s_df, s_s2n = felis.load_spectrum(spectrum,verbose=verbose)


    if waverange is not None:
        if verbose: print(' - Limiting the wavelength range to cross-correlate over ')
        goodent = np.where( (s_wave > waverange[0]) & (s_wave < waverange[1]) )
        s_wave, s_flux, s_df, s_s2n = s_wave[goodent], s_flux[goodent], s_df[goodent], s_s2n[goodent]

    if verbose: print(' - Estimating continuum from spectrum (median value) in observed units')
    if spec_median_sub:
        continuumval = np.median(s_flux)
    else:
        continuumval = 0.0

    if z_restframe is not None:
        if verbose: print(' - Shifting spectrum to rest-frame using the redshift '+str(z_restframe))
        s_wave, s_flux, s_df, s_s2n = s_wave / (1+z_restframe), s_flux * (1+z_restframe), s_df * (1+z_restframe), s_s2n
    else:
        z_restframe=0.0

    if verbose: print(' - Removing continuum from spectrum ')
    s_flux       = s_flux - continuumval * (1.0+z_restframe)

    if verbose: print(' - Interpolate template to spectrums wavelength resolution before cross-correlating')
    t_wave_init, t_flux_init, t_df_init, t_s2n_init = felis.load_spectrum(template,verbose=verbose)
    func       = scipy.interpolate.interp1d(t_wave_init,t_flux_init,kind='linear',fill_value="extrapolate")
    t_flux     = func(s_wave)
    t_flux[t_flux < min_template_level] = 0.0

    if verbose: print(' - Normalizing total template flux to 1')
    if len(t_flux[t_flux != 0]) == 0:
        if verbose: print('   FELIS WARNING All interpolated template pixels are 0.0')
    else:
        temp_int = np.trapz(t_flux,s_wave)
        t_flux   = t_flux / temp_int

    Npix = len(s_flux)

    template_triplelength = np.zeros(3*Npix)
    template_triplelength[0:Npix] = t_flux

    ccresults = np.zeros([Npix,5])

    for ii in np.arange(Npix):
        rollsize       = int(ii+np.floor(Npix/2.))
        template_shift = np.roll(template_triplelength, rollsize)[Npix:-Npix]

        normlim   = 1e-10
        temp_int  = np.trapz(template_shift,s_wave)
        sumdiff   = np.abs(temp_int-1.0)
        if sumdiff > normlim: # make sure "half-template" shifts are also normalized
            if temp_int == 0: # only attempt normalization if template is not just zeros
                pass
            else:
                template_shift   = template_shift / temp_int

        try:
            flux_scale, flux_scale_variance, S2N, chi2_min, NgoodentChi2 = \
                felis.minimize_chi2(s_flux,s_df,template_shift,temp_wave=s_wave,verbose=False)

            ccresults[ii,:] = flux_scale, flux_scale_variance, S2N, chi2_min, NgoodentChi2
        except:
            print(' ERROR: Problems in minimizing Chi**2 with felis.minimize_chi2() while cross-correlating. '
                  'Stopping for further investigation')

            pdb.set_trace()

    max_S2N     = np.max(ccresults[:,2])
    max_S2N_ent = np.where(ccresults[:,2] == max_S2N)[0][0]
    max_wave    = s_wave[max_S2N_ent]

    dlam     = np.median(np.diff(s_wave))
    t_wave   = np.arange(np.min(t_wave_init),np.max(t_wave_init),dlam)
    t_cenent = int(np.ceil(len(t_wave)/2))
    max_z    = (max_wave*(z_restframe+1.0) / t_wave[t_cenent]) - 1.0

    return  s_wave, ccresults, max_S2N, max_z, continuumval

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
def save_spectrum(outfile,wave,flux,fluxerr,headerinfo=None,waveunits='ANGSTROMS',fluxunits='',overwrite=False,verbose=True):
    """
    Function saving a spectrum (or template) to a fits file

    --- INPUT ---
    outfile         Name of output file to generate
    wave            Wavelength in units of Angstrom
    flux            Flux to store to fits file
    fluxerr         Error un flux
    headerinfo      To add info to the header provide it to this keyword as a dictionary on the format:
                       headerinfo = {'KEYNAME1':[VALUE1,INFOCOMMENT1], 'KEYNAME2':[VALUE2,INFOCOMMENT2]}
    waveunits       Units of wave vector to store in fits header
    fluxunits       Units of flux to store in fits header
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
    if verbose: print(' - Saving wavelength and flux values to \n   '+outfile)
    mainHDU = fits.PrimaryHDU()       # primary HDU
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    c1 = fits.Column(name='wave',      format='D', unit=waveunits, array=wave)
    c2 = fits.Column(name='flux',      format='D', unit=fluxunits, array=flux)
    c3 = fits.Column(name='fluxerror', format='D', unit=fluxunits, array=fluxerr)
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
    pversion = sys.version_info[0]
    with open(picklefile, 'rb') as f:
        if pversion == 2:
            return pickle.load(f) # for python 2.X
        elif pversion == 3:
            return pickle.load(f, encoding='latin1') # for python 3.X
        else:
            sys.exit('   felis.load_picklefile() saw unknown version of python: version = '+str(pversion))
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def getresult4maxS2N(dictionary,dictionarykey,zspecISzLya=False):
    """
    Function returning results at maxmimum S/N for the template with maximum S/N overall, for a given
    dictionary keyword of the pickled FELIS output.

    --- INPUT ---
    dictionary      = Dictionary containing FELIS outputs, i.e. the output from felis.load_picklefile()
    dictionarykey   = The object/spectrum (dictionary key) to return results for
    zspecISzLya         If picklefile contains zLya keyword insead of zspec (was generated before 180912), set this
                        keyword to True to enable proper handling of the dictionary keywords.

    """
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Ensure compatibility with FELIS output dictionaries from before 180912
    zkey = 'zspec'
    if zspecISzLya: zkey = 'zLya'
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    keys = dictionary.keys()
    if dictionarykey in keys:
        keydic        = dictionary[dictionarykey]
        S2Nmaxvec     = keydic['S2NCCmaxvec']
        S2Nmax        = np.max(S2Nmaxvec)
        templateent   = np.where(S2Nmaxvec == S2Nmax)[0]
        if len(templateent) > 1:
            print('    FELIS WARNING felis.getresult4maxS2N(): \n    Multiple templates ('+str(len(templateent))+' of '+
                  str(len(S2Nmaxvec))+' templates) share the same maximum S/N ('+str(S2Nmax)+
                  ') of cross-correlation for \n   '+dictionarykey)
            print('   Will return the results for the first template')
            templateent = np.array([templateent[0]])

        template      = keydic['templatevec'][templateent][0]
        zS2Nmax       = keydic['zCCmaxvec'][templateent][0]
        zspec         = keydic[zkey]
        vshift_match  = 299792.458 * (zspec - zS2Nmax)/(1.0+zS2Nmax) # cf. Erb+2014

        if 'vshift' in keydic.keys():
            vshift_intr   = keydic['vshift'] # velocity shift stored in dictionary with felis.match_templates2specs()
        else:
            vshift_intr   = -99

        S2Nvec        = keydic['ccresultsarr_S2N'][templateent][0]
        S2Nvecmax     = np.where(S2Nvec == np.max(S2Nvec))[0]

        flux          = keydic['ccresultsarray_flux'][templateent,S2Nvecmax][0]
        fluxerr       = np.sqrt(keydic['ccresultsarray_variance'][templateent,S2Nvecmax][0])
        Ngoodent      = keydic['ccresultsarr_Ngoodent'][templateent,S2Nvecmax][0]
        chi2          = keydic['ccresultsarr_chi2'][templateent,S2Nvecmax][0]
    else:
        template      = 'None'
        vshift_intr   = -99.0
        vshift_match  = -99.0
        S2Nmax        = -99.0
        flux          = -99.0
        fluxerr       = -99.0
        Ngoodent      = -99.0
        chi2          = -99.0

    return template, vshift_intr, vshift_match, flux, fluxerr, S2Nmax, Ngoodent, chi2, zspec, zS2Nmax
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
