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
    '''
    
    log.info('### Begin main ! ')

    phot0, err_phot0, sptype, tab0 = get_photometry.main(name)
    mag0 = [phot0['J'], phot0['H'], phot0['K']]
    
    wave, F_nu = read_stellar_template.main(sptype, library=library)

    model_mag0 = spec_convolve.main(wave, F_nu, AB=False)

    mag_diff     = mag0 - model_mag0
    avg_mag_diff = np.average(mag_diff)

    F_nu /= 10**(avg_mag_diff/2.5)

    log.info('### End main ! ')
#enddef

