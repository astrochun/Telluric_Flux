"""
get_photometry
==============

Get photometry from Simbad for telluric stars
"""

import matplotlib.pyplot as plt

from astropy import log

from astroquery.simbad import Simbad

add_fields = ['sptype', 'flux(J)', 'flux_error(J)', 'flux(H)', 'flux_error(H)',
              'flux(K)', 'flux_error(K)']
for field0 in add_fields:
    if not any(field0 in sfield for sfield in Simbad.get_votable_fields()):
        Simbad.add_votable_fields(field0)

def main(name):

    '''
    Main function to query Simbad to get photometry for telluric stars

    Parameters
    ----------

    Returns
    -------

    Notes
    -----
    Created by Chun Ly, 17 April 2018
    '''

    log.info('### Begin main ! ')

    tab0 = Simbad.query_object(name)

    phot0     = {'J': tab0['FLUX_J'][0], 'H': tab0['FLUX_H'][0],
                 'K': tab0['FLUX_K'][0]}
    err_phot0 = {'J': tab0['FLUX_ERROR_J'][0], 'H': tab0['FLUX_ERROR_H'][0],
                 'K': tab0['FLUX_ERROR_K'][0]}
    sptype    = tab0['SP_TYPE'][0].lowercase()

    log.info('### End main ! ')

    return phot0, err_phot0, sptype, tab0
#enddef

