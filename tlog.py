"""
tlog
====

Routine to log information to stdout and ASCII file

To execute:
from Telluric_Flux import tlog

logfile = '/path/to/telluric/data/telluric_flux.log'
mylogger = tlog.log0(logfile)._get_logger()
"""

import logging, sys

formatter = logging.Formatter('%(asctime)s - %(module)12s.%(funcName)20s - %(levelname)s: %(message)s')

# set up logging to STDOUT for all levels INFO and higher
sh = logging.StreamHandler(sys.stdout)
sh.setLevel(logging.INFO)
sh.setFormatter(formatter)

class log0:
    '''
    Main class to log information to stdout and ASCII file

    To execute:
    logfile = '/path/to/telluric/data/telluric_flux.log'
    mylogger = log0(logfile)._get_logger()

    Parameters
    ----------

    Returns
    -------

    Notes
    -----
    Created by Chun Ly, 18 April 2018
    '''

    def __init__(self,file):
        self.LOG_FILENAME = file
        self._log = self._get_logger()

    def _get_logger(self):
        loglevel = logging.INFO
        log = logging.getLogger(self.LOG_FILENAME) # + Mod on 14/12/2017
        if not getattr(log, 'handler_set', None):
            log.setLevel(logging.INFO)
            sh = logging.StreamHandler()
            sh.setFormatter(formatter)
            log.addHandler(sh)

            fh = logging.FileHandler(self.LOG_FILENAME)
            fh.setLevel(logging.INFO)
            fh.setFormatter(formatter)
            log.addHandler(fh)

            log.setLevel(loglevel)
            log.handler_set = True
        return log
#enddef
