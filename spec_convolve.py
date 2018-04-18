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

co_dir0   = dirname(__file__)
filt_dir0 = co_dir0 + '/filters/'

def main(wave, F_nu):

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

    Notes
    -----
    Created by Chun Ly, 18 April 2018
    '''

    log.info('### Begin main ! ')

    # Read in 2MASS J filter
    filt_files = glob(filt_dir0+'2MASS_?.txt')

    for ff in range(len(filt_files)):
        filt_tab = asc.read(filt_files[ff])
        filt_wave = filt_tab['col1'] * 1E4

        f = interp1d(filt_wave, filt_tab['col2'])
        wave_interp = f(wave)

    log.info('### End main ! ')
#enddef

