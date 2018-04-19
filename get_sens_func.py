"""
get_sens_func
=============

Get sensitivity function (erg/s/cm2/AA -> DN/s) using telluric star spectra
"""

from astropy.io import ascii as asc
from astropy.io import fits

import numpy as np

import matplotlib.pyplot as plt
import glob

from astropy import log

from . import read_stellar_template, spec_convolve, get_photometry

c_Ang = spec_convolve.c_Ang

c_Jy = 1e-23 # erg/s/cm2/Hz to Janskies

Ang_micron = 1e4 # Angstrom to micron conversion

def main(name, library='Pickles'):

    '''
    Main function to call all functions to get stellar spectra template,
    normalize by photometry, and then determine sensitivity function based
    on 1-D spectra

    Parameters
    ----------

    Returns
    -------

    Notes
    -----
    Created by Chun Ly, 18 April 2018
     - Get filter center and width from spec_convolve()
     - Plot photometry-normalized F_nu vs lambda
     - Call spec_convolve.get_vega_fluxes()
     - Plot photometric data on photometry-normalized F_nu
     - Plot aesthetics: Jy for F_nu, errorbars, output to file
     - Plot aesthetics: Plot F_lam in bottom panel; Angstroms -> microns
    '''
    
    log.info('### Begin main ! ')

    phot0, err_phot0, sptype, tab0 = get_photometry.main(name)
    mag0 = np.array([phot0['J'], phot0['H'], phot0['K']])

    # Get Vega fluxes to convert from Vega mag to F_nu | + on 19/04/2018
    Vega_fluxes = spec_convolve.get_vega_fluxes()
    flux0       = 10**(-0.4*mag0) * Vega_fluxes

    wave, F_nu = read_stellar_template.main(sptype, library=library)

    model_mag0, wave_cen0, FWHM0 = spec_convolve.main(wave, F_nu, AB=False)
    print wave_cen0, FWHM0

    mag_diff     = mag0 - model_mag0
    avg_mag_diff = np.average(mag_diff)

    F_nu /= 10**(avg_mag_diff/2.5)

    F_lam = F_nu * c_Ang / (wave**2)

    fig, ax = plt.subplots(nrows=2) # + on 19/04/2018

    # + on 19/04/2018
    F_Jy = F_nu/c_Jy
    ax[0].semilogy(wave/Ang_micron, F_Jy, 'b-', zorder=1)

    # Plot photometric data on photometry-normalize F_nu | + on 19/04/2018
    ax[0].errorbar(wave_cen0/Ang_micron, flux0/c_Jy, xerr=FWHM0/2.0/Ang_micron,
                   ecolor='black', elinewidth=1.5, capsize=3.0, capthick=1.5,
                   zorder=2, fmt='none')
    ax[0].scatter(wave_cen0/Ang_micron, flux0/c_Jy, 50, marker='o', color='red',
                  edgecolor='black', zorder=2)
    ax[0].set_xlabel('')
    ax[0].set_xticklabels([])
    ax[0].set_ylabel(r'$F_{\nu}$ [Jy]')
    ax[0].minorticks_on()

    # + on 19/04/2018
    xlim = np.array([min(wave)-100,max(wave)+100])/Ang_micron
    ax[0].set_xlim(xlim)

    # + on 19/04/2018
    non_zero = np.where(F_Jy != 0)[0]
    ylim = [10**(np.floor(np.log10(min(F_Jy[non_zero])))), 2*max(F_Jy)]
    ax[0].set_ylim(ylim)

    # Plot F_lam | + on 19/04/2018
    ax[1].semilogy(wave/Ang_micron, F_lam, 'b-')
    ax[1].set_xlim(xlim)

    non_zero = np.where(F_lam != 0)[0]
    ylim = [10**(np.floor(np.log10(min(F_lam[non_zero])))), 2*max(F_lam)]
    ax[1].set_ylim(ylim)

    ax[1].set_xlabel(r'Wavelengths [$\mu$m]')
    ax[1].set_ylabel(r'$F_{\lambda}$ [erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$]')
    ax[1].minorticks_on()

    # + on 19/04/2018
    plt.subplots_adjust(left=0.1, right=0.98, bottom=0.08, top=0.99,
                        hspace=0.03)
    fig.savefig(name+'_spec_model.pdf')

    log.info('### End main ! ')
#enddef

