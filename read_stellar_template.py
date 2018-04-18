"""
read_stellar_template
====

Provide description for code here.
"""

from os.path import dirname

from os.path import exists
from astropy.io import ascii as asc
from astropy.io import fits

import numpy as np

import matplotlib.pyplot as plt
from glob import glob

from astropy import log
import astropy.constants as const
import astropy.units as u

co_dir0   = dirname(__file__)
temp_dir0 = co_dir0 + '/templates/'

def main(stellar_type, library='Pickles'):

    '''
    Main function to read in template

    Parameters
    ----------
    stellar_type : str
      String for type. Lower case (e.g., a0v)

    library : str
      Stellar library to use. Options are: 'Pickles'
      Default: Pickles

    silent : boolean
      Turns off stdout messages. Default: False

    verbose : boolean
      Turns on additional stdout messages. Default: True

    Returns
    -------
    wave : np.array
      Array of wavelengths in Angstroms

    F_nu : np.array
      Array of F_nu in erg/s/cm2/Hz (not normalized to photometry)

    Notes
    -----
    Created by Chun Ly, 17 April 2018
     - Return wave and F_nu
     - Bug fix: CDELT -> CDELT1
    '''

    log.info('### Begin main !')

    temp_dir = temp_dir0 + library + '/'
    temp_file = glob(temp_dir + '*'+stellar_type+'.fits')
    if len(temp_file) == 1: temp_file = temp_file[0]
    log.info('### '+temp_file)

    log.info('### Reading : '+temp_file)
    F_lam, hdr = fits.getdata(temp_file, header=True)

    crval1 = hdr['CRVAL1']
    cdelt1 = hdr['CDELT1']
    NX     = hdr['NAXIS1']
    wave = crval1 + cdelt1 * np.arange(NX)

    # erg/s/cm2/Ang -> erg/s/cm2/Hz
    F_nu = F_lam * (wave**2) / const.c.to(u.Angstrom/u.s).value

    return wave, F_nu
    
    log.info('### End main !')
#enddef

