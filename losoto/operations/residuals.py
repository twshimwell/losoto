#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Residual operation for LoSoTo

# This operation subtract a clock and/or tec from a phase

import logging
from losoto.operations_lib import *

logging.debug('Loading RESIDUAL module.')

def run( step, parset, H ):
    """
    subtract a clock and/or tec from a phase.
    """
    import numpy as np
    from losoto.h5parm import solFetcher, solWriter

    soltabs = getParSoltabs( step, parset, H )
    soltabsToSub = parset.getStringVector('.'.join(["LoSoTo.Steps", step, "Sub"]), [] )

    # select only clock/tec tables for subtraction
#    sfss = []
#    for soltabToSub in soltabsToSub
#        sf = solFetcher(soltabToSub)
#        if sf.getType() != 'clock' or sf.getType() != 'tec':
#            logging.warning(soltab._v_name+' is not of type clock/tec and cannot be subtracted, ignore.')
#        else:
#            sfss.append(sf)

    for soltab in openSoltabs( H, soltabs ):
        logging.info("--> Working on soltab: "+soltab._v_name)

        sf = solFetcher(soltab)
        sw = solWriter(soltab)
        if sf.getType() != 'phase':
            logging.warning(soltab._v_name+' is not of type phase, ignore.')
            continue

        # the only return axes is freq, slower but better code
        for vals, weights, coord, selection in sf.getValuesIter(returnAxes='freq', weight = True):

            for soltabToSub in soltabsToSub:
                sfs = solFetcher(soltabToSub)
                # restrict to phase coordinates
                newCoord = {}
                for axisName in coord.keys():
                    if axisName in sfs.getAxesNames(): newCoord[axisName] = coord[axisName]
                sfs.setSelection(**newCoord)
                valsSub = np.squeeze(sfs.getValues(retAxesVals=False, weight=False))
                weightsSub = np.squeeze(sfs.getValues(retAxesVals=False, weight=True))

                if sfs.getTyep == 'clock':
                    logging.info('Subtracting CLOCK table: '+soltabToSub._v_name)
                    vals -= 2. * np.pi * valsSub * coord['freq']

                elif sfs.getTyep == 'tec':
                    logging.info('Subtracting TEC table: '+soltabToSub._v_name)
                    vals -= -8.44797245e9 * valsSub / coord['freq']

                else:
                    logging.warning(soltabToSub._v_name+' is not of type clock/tec and cannot be subtracted, ignore.')

            # Write back an AND of the flags
            sw.selection = selection
            sw.setValues(vals)
        sw.addHistory('SUBTRACTED with tables '+' '.join(soltabsToSub))
        
    return 0
