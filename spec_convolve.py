"""
spec_convolve
=============

Convolve filter bandpasses with spectra
"""

from os.path import dirname

import numpy as np
from scipy.interpolate import interp1d

from glob import glob

import matplotlib.pyplot as plt

import astropy.io.ascii as asc
from astropy import log
import astropy.constants as const
import astropy.units as u

co_dir0   = dirname(__file__)
filt_dir0 = co_dir0 + '/filters/'
temp_dir0 = co_dir0 + '/templates/'

c_Ang = const.c.to(u.Angstrom/u.s).value

# Read in Vega SED  | Moved up on 18/04/2018
Vega_file = temp_dir0+'Vega.sed'
Vega_tab  = asc.read(Vega_file)

def get_filter_center(wave, filt_profile):

    dx     = wave[1]-wave[0]
    top    = np.sum(filt_profile * wave * dx)
    bottom = np.sum(filt_profile * dx)
    cen_wave = top/bottom

    top    = np.sum(filt_profile * (wave-cen_wave)**2 * dx)
    bottom = np.sum(filt_profile * dx)
    FWHM_wave = 2 * np.sqrt(2*np.log(2)) * np.sqrt(top/bottom)
    return cen_wave, FWHM_wave
#enddef

def main(wave, F_nu, AB=False):

    '''
    Main function to convolve spectral templates with filters

    Parameters
    ----------

    silent : boolean
      Turns off stdout messages. Default: False

    verbose : boolean
      Turns on additional stdout messages. Default: True

    Returns
    -------
    mag0_arr : np.array
      Array of magnitudes

    Notes
    -----
    Created by Chun Ly, 17 April 2018

    Modified by Chun Ly, 18 April 2018
     - Determine fluxes and AB magnitudes from convolution with filters
     - Determine fluxes and Vega magnitudes from convolution with filters
     - Change waveband order for 2MASS filters
     - Call get_filter_center()
     - Moved Vega SED definition outside main()
    '''

    log.info('### Begin main ! ')

    # Read in 2MASS J filter
    filt_files = [filt_dir0+'2MASS_'+bands+'.txt' for bands in ['J','H','K']]

    nu = c_Ang/wave

    mag0_arr  = np.zeros(len(filt_files))
    cen0_arr  = np.zeros(len(filt_files))
    FWHM0_arr = np.zeros(len(filt_files))

    for ff in range(len(filt_files)):
        filt_tab = asc.read(filt_files[ff])
        log.info('## Reading : '+filt_files[ff])
        filt_wave = filt_tab['col1'] * 1E4

        cen0_arr[ff], FWHM0_arr[ff] = get_filter_center(filt_wave, filt_tab['col2'])

        f = interp1d(filt_wave, filt_tab['col2'], bounds_error=False,
                     fill_value = 0.0)
        filt_interp = f(wave)

        d_nu       = np.zeros(len(wave))
        d_nu[0:-1] = nu[1:] - nu[0:-1]
        d_nu[-1]   = d_nu[-2]

        top    = np.sum(F_nu * filt_interp * d_nu)
        if AB == True:
            bottom = np.sum(filt_interp * d_nu)
            mag0_arr[ff] = -2.5*np.log10(top/bottom) - 48.6
        else: # + on 18/04/2018
            V_wave = Vega_tab['col1']
            V_Flam = Vega_tab['col2']
            V_Fnu  = V_Flam * (V_wave**2) / c_Ang
            f_Vega = interp1d(V_wave, V_Fnu)
            V_Fnu_interp = f_Vega(wave)

            bottom = np.sum(filt_interp * V_Fnu_interp * d_nu)
            mag0_arr[ff] = -2.5*np.log10(top/bottom)
    #endfor

    log.info('### End main ! ')
    return mag0_arr, cen0_arr, FWHM0_arr

#enddef

def get_vega_fluxes():
    '''
    Compute Vega F_nu fluxes in each waveband. This is to normalize magnitudes
    to spectra

    Parameters
    ----------
    None.

    Returns
    -------
    flux0_arr : np.array
      Array of fluxes

    Notes
    -----
    Created by Chun Ly, 19 April 2018
    '''

    # Read in 2MASS J filter
    filt_files = [filt_dir0+'2MASS_'+bands+'.txt' for bands in ['J','H','K']]

    flux0_arr = np.zeros(len(filt_files))

    for ff in range(len(filt_files)):
        filt_tab = asc.read(filt_files[ff])
        log.info('## Reading : '+filt_files[ff])
        filt_wave = filt_tab['col1'] * 1E4
        nu = c_Ang/filt_wave

        V_wave = Vega_tab['col1']
        V_Flam = Vega_tab['col2']
        V_Fnu  = V_Flam * (V_wave**2) / c_Ang
        f_Vega = interp1d(V_wave, V_Fnu)
        V_Fnu_interp = f_Vega(filt_wave)

        d_nu       = np.zeros(len(filt_wave))
        d_nu[0:-1] = nu[1:] - nu[0:-1]
        d_nu[-1]   = d_nu[-2]

        top    = np.sum(V_Fnu_interp * filt_tab['col2'] * d_nu)
        bottom = np.sum(filt_tab['col2'] * d_nu)
        flux0_arr[ff] = top/bottom
    #endfor
    return flux0_arr
#enddef
