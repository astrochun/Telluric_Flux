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

c_Ang = const.c.to(u.Angstrom/u.s).value

def main(wave, F_nu, AB=False):

    '''
    Provide explanation for function here.

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
    '''

    log.info('### Begin main ! ')

    # Read in 2MASS J filter
    filt_files = glob(filt_dir0+'2MASS_?.txt')

    nu = c_Ang/wave

    mag0_arr = np.zeros(len(filt_files))

    for ff in range(len(filt_files)):
        filt_tab = asc.read(filt_files[ff])
        filt_wave = filt_tab['col1'] * 1E4

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

    #endfor

    log.info('### End main ! ')
    return mag0_arr

#enddef

