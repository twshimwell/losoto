#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This tool is used to convert an H5parm file to parmdb format by writing to
# existing parmdb instrument table(s).
#
# It handles Gain/DirectionalGain/RotationAngle/CommonRotationAngle/CommonScalarPhase solution types.
_author = "Francesco de Gasperin (fdg@hs.uni-hamburg.de), David Rafferty (drafferty@hs.uni-hamburg.de)"

import sys, os, glob, re, time
import numpy as np
import shutil
import logging
import pyrap.tables as pt
import lofar.parmdb
from losoto import _version
from losoto import _logging
from losoto.h5parm import h5parm, solWriter, solFetcher
try:
    import progressbar
except ImportError:
    import losoto.progressbar as progressbar


def getSoltypeFromSolTabs(solTabs):
    """Return a list of parmdb solution type given for input tables

    This function is basically the reverse of parmdbToAxes() in
    H5parm_importer.py.

    solTabs - solution tables returned by h5parm.getSoltabs()

    """
    solTabList = []

    for name, st in solTabs.iteritems():
        if st._v_title == 'amplitude':
            if dir == 'pointing':
                solType == 'Gain'
            else:
                solType == 'DirectionalGain'
        elif st._v_title == 'phase':
            if dir == 'pointing':
                solType == 'Gain'
            else:
                solType == 'DirectionalGain'
        elif st._v_title == 'rotation':
            if dir == 'pointing':
                solType == 'CommonRotationAngle'
            else:
                solType == 'RotationAngle'
        elif st._v_title == 'scalarphase':
            if dir == 'pointing':
                solType == 'CommonScalarPhase'
            else:
                solType == 'ScalarPhase'
        elif st._v_title == 'scalaramplitude':
            if dir == 'pointing':
                solType == 'CommonScalarAmplitude'
            else:
                solType == 'ScalarAmplitude'
        elif st._v_title == 'clock':
            solType == 'Clock'
        elif st._v_title == 'tec':
            solType == 'TEC'
        elif st._v_title == 'rotationmeasure':
            solType == 'RotationMeasure'

    return solTabList


def makeTECparmdb(H, solset, TECsolTab, timewidths, freq, freqwidth):
    """Returns TEC screen parmdb parameters

    H - H5parm object
    solset - solution set with TEC screen parameters
    TECsolTab = solution table with tecscreen values
    timewidths - time widths of output parmdb
    freq - frequency of output parmdb
    freqwidth - frequency width of output parmdb
    """
    global ipbar, pbar

    station_dict = H.getAnt(solset)
    station_names = station_dict.keys()
    station_positions = station_dict.values()
    source_dict = H.getSou(solset)
    source_names = source_dict.keys()
    source_positions = source_dict.values()

    tec_sf = solFetcher(TECsolTab)
    tec_screen, axis_vals = tec_sf.getValues()
    times = axis_vals['time']
    beta = TECsolTab._v_attrs['beta']
    r_0 = TECsolTab._v_attrs['r_0']
    height = TECsolTab._v_attrs['height']
    order = TECsolTab._v_attrs['order']
    pp = tec_sf.t.piercepoint

    N_sources = len(source_names)
    N_times = len(times)
    N_freqs = 1
    N_stations = len(station_names)
    N_piercepoints = N_sources * N_stations

    freqs = freq
    freqwidths = freqwidth
    parms = {}
    v = {}
    v['times'] = times
    v['timewidths'] = timewidths
    v['freqs'] = freqs
    v['freqwidths'] = freqwidths

    for station_name in station_names:
        for source_name in source_names:

            v['values'] = np.zeros((N_times, N_freqs), dtype=np.double)
            parmname = 'Piercepoint:X:%s:%s' % (station_name, source_name)
            parms[parmname] = v.copy()

            v['values'] = np.zeros((N_times, N_freqs), dtype=np.double)
            parmname = 'Piercepoint:Y:%s:%s' % (station_name, source_name)
            parms[parmname] = v.copy()

            v['values'] = np.zeros((N_times, N_freqs), dtype=np.double)
            parmname = 'Piercepoint:Z:%s:%s' % (station_name, source_name)
            parms[parmname] = v.copy()

            v['values'] = np.zeros((N_times, N_freqs), dtype=np.double)
            parmname = 'TECfit_white:%s:%s' % (station_name, source_name)
            parms[parmname] = v.copy()

            v['values'] = np.zeros((N_times, N_freqs), dtype=np.double)
            parmname = 'TECfit_white:0:%s:%s' % (station_name, source_name)
            parms[parmname] = v.copy()

            v['values'] = np.zeros((N_times, N_freqs), dtype=np.double)
            parmname = 'TECfit_white:1:%s:%s' % (station_name, source_name)
            parms[parmname] = v.copy()

    for k in range(N_times):
        D = np.resize(pp[k, :, :], (N_piercepoints, N_piercepoints, 3))
        D = np.transpose(D, ( 1, 0, 2 )) - D
        D2 = np.sum(D**2, axis=2)
        C = -(D2 / (r_0**2))**(beta / 2.0) / 2.0
        tec_fit_white = np.dot(np.linalg.inv(C),
            tec_screen[:, k, :].reshape(N_piercepoints))
        pp_idx = 0
        for src, source_name in enumerate(source_names):
            for sta, station_name in enumerate(station_names):

                parmname = 'Piercepoint:X:%s:%s' % (station_name, source_name)
                parms[parmname]['values'][k, 0] = pp[k, pp_idx, 0]

                parmname = 'Piercepoint:Y:%s:%s' % (station_name, source_name)
                parms[parmname]['values'][k, 0] = pp[k, pp_idx, 1]

                parmname = 'Piercepoint:Z:%s:%s' % (station_name, source_name)
                parms[parmname]['values'][k, 0] = pp[k, pp_idx, 2]

                parmname = 'TECfit_white:%s:%s' % (station_name, source_name)
                parms[parmname]['values'][k, 0] = tec_fit_white[pp_idx]

                parmname = 'TECfit_white:0:%s:%s' % (station_name, source_name)
                parms[parmname]['values'][k, 0] = tec_fit_white[pp_idx]

                parmname = 'TECfit_white:1:%s:%s' % (station_name, source_name)
                parms[parmname]['values'][k, 0] = tec_fit_white[pp_idx]

                pp_idx += 1
        pbar.update(ipbar)
        ipbar += 1

    time_start = times[0] - timewidths[0]/2
    time_end = times[-1] + timewidths[-1]/2

    v['times'] = np.array([(time_start + time_end) / 2])
    v['timewidths'] = np.array([time_end - time_start])

    v_r0 = v.copy()
    v_r0['values'] = np.array(r_0, dtype=np.double, ndmin=2)
    parms['r_0'] = v_r0

    v_beta = v.copy()
    v_beta['values'] = np.array(beta, dtype=np.double, ndmin=2)
    parms['beta'] = v_beta

    v_height = v.copy()
    v_height['values'] = np.array(height, dtype=np.double, ndmin=2)
    parms['height'] = v_height

    return parms


if __name__=='__main__':
    # Options
    import optparse
    opt = optparse.OptionParser(usage='%prog <H5parm filename> <output parmdb>\n'+
        _author, version='%prog '+_version.__version__)
    opt.add_option('-v', '--verbose', help='Go VeRbOsE!',
        action='store_true', default=False)
    opt.add_option('-s', '--solset', help='Name of solution set to export '
        '(default=sol000)', type='string', default='sol000')
    opt.add_option('-t', '--soltab', help='Solution tables to export; e.g., '
        '"amplitude000, phase000" (default=all)', type='string', default='all')
    opt.add_option('-c', '--clobber', help='Clobber exising files '
        '(default=False)', action='store_true', default=False)
    (options, args) = opt.parse_args()
    global ipbar, pbar

    # Check options
    if len(args) != 2:
        opt.print_help()
        sys.exit()
    if options.verbose: _logging.setLevel("debug")

    # Check input H5parm file
    h5parmFile = args[0]
    if not os.path.exists(h5parmFile):
        logging.critical('Input H5parm file not found.')
        sys.exit(1)
    logging.info("Input H5parm filename = "+h5parmFile)

    # Open the H5parm file and get solution set names
    h5parm_in = h5parm(h5parmFile, readonly = True)
    solsetNames = h5parm_in.getSolsets()

    # Check input solution set name
    solsetName = options.solset
    if solsetName not in solsetNames:
        logging.critical('The solution set "'+solsetName+'" was not found in input H5parm file.')
        sys.exit(1)
    logging.info("Solution set name = "+solsetName)
    solset = h5parm_in.getSolset(solsetName)

    # Check output parmdb file
    globaldbFile = args[1]
    logging.info("Output globaldb/SB filename = "+globaldbFile)
    if not os.path.exists(globaldbFile):
        out_globaldbFile_exists = False
        logging.info('Output globaldb/SB file not found. It will be created.')
    else:
        out_globaldbFile_exists = True

    # Find solution table types
    solTabs = h5parm_in.getSoltabs(solset)
    if options.soltab != 'all':
        soltabs_to_use = [s.strip() for s in options.soltab.split(',')]
        logging.info('Using solution tables: {0}'.format(soltabs_to_use))
        solTabs_filt = {}
        for s, v in solTabs.iteritems():
            if s in soltabs_to_use:
                solTabs_filt[s] = v
        for s in soltabs_to_use:
            if s not in solTabs_filt.keys():
                logging.warning('Solution table {0} not found in input H5parm file.'.format(s))
        solTabs = solTabs_filt
    if len(solTabs) == 0:
        logging.critical('No solution tables found matching input criteria.')
        sys.exit(1)
    solTypes = getSoltypeFromSolTabs(solTabs)
    if len(solTypes) is None:
        logging.warning("No valid solutions found in solution set {0}.".format(solsetName))
        sys.exit(1)

    solTypes = list(set(solTypes))
    logging.info('Exporting the following solution types: '+', '.join(solTypes))

    # For each solType, select appropriate solutions and construct
    # the dictionary to pass to pdb.addValues()
    len_sol = {}
    for solType in solTypes:
        if solType != 'TECScreen':
            len_sol[solType] = len(pdb.getNames(solType+':*'))
        else:
            tec_sf = solFetcher(st_tec)
            N_times = tec_sf.getAxisLen(axis='time')
            len_sol[solType] = N_times

    cachedSolTabs = {}
    for instrumentdbFile in instrumentdbFiles:
        out_instrumentdbFile = out_globaldbFile + '/' + outroot + '_' + instrumentdbFile.split('/')[-1]
        logging.info('Filling '+out_instrumentdbFile+':')

        # Remove existing instrumentdb (if clobber) and create new one
        if os.path.exists(out_instrumentdbFile):
            if options.clobber:
                shutil.rmtree(out_instrumentdbFile)
            else:
                logging.critical('Output instrumentdb file exists and '
                    'clobber = False.')
                sys.exit(1)
        pdb_out = lofar.parmdb.parmdb(out_instrumentdbFile+'/', create=True)

        pbar = progressbar.ProgressBar(maxval=sum(len_sol.values())).start()
        ipbar = 0

        pdb_in = lofar.parmdb.parmdb(instrumentdbFile)

        # Add default values and steps
        DefValues = pdb_in.getDefValues()
        for k, v in DefValues.iteritems():
            pdb_out.addDefValues({k: pdb.makeDefValue(v.item(0))})
        pdb_out.setDefaultSteps(pdb_in.getDefaultSteps())

        for solType in solTypes:
            if len_sol[solType] == 0: continue

            if solType != 'TECScreen':
                solEntries = pdb_in.getNames(solType+':*')
                data = pdb_in.getValuesGrid(solType+':*')
                data_out = data.copy()
                for solEntry in solEntries:

                    pol, dir, ant, parm = parmdbToAxes(solEntry)
                    solTabList = getSoltabFromSolType(solType, solTabs, parm=parm)
                    if solTabList is None:
                        continue
                    if len(solTabList) > 1:
                        logging.warning('More than one solution table found in H5parm '
                                'matching parmdb entry "'+solType+'". Taking the first match: '+str(solTabList[0])+'.')
                    solTab = solTabList[0]

                    # search in the cache for open soltab
                    if not solTab._v_title in cachedSolTabs:
                        sf = solFetcher(solTab, useCache=True)
                        cachedSolTabs[solTab._v_title] = sf
                    else:
                        sf = cachedSolTabs[solTab._v_title]

                    freqs = data[solEntry]['freqs']
                    times = data[solEntry]['times']
                    parms = {}
                    if 'ant' in sf.getAxesNames(): parms['ant'] = [ant]
                    if 'pol' in sf.getAxesNames(): parms['pol'] = [pol]
                    if 'dir' in sf.getAxesNames(): parms['dir'] = [dir]
                    if 'freq' in sf.getAxesNames(): parms['freq'] = freqs.tolist()
                    # workaround for bbs and ndppp dealing differently with the last time slot when #timeslots%ntime != 0
                    # NDPPP has all intervals the same
                    # BBS has a maller interval in the last timeslot which is compensated here
                    if times[-1] - times[-2] < times[-2] - times[-3]: times[-1] = times[-2] + (times[-2] - times[-3])
                    if 'time' in sf.getAxesNames(): parms['time'] = {'min':np.min(times-0.1), 'max':np.max(times+0.1)}
                    sf.setSelection(**parms)

                    # If needed, convert Amp and Phase to Real and Imag
                    if parm == 'Real':
                        solTabList = getSoltabFromSolType(solType, solTabs, parm='phase')
                        if not solTabList[0]._v_title in cachedSolTabs:
                            sf_phase = solFetcher(solTabList[0], useCache=True)
                            cachedSolTabs[solTabList[0]._v_title] = sf_phase
                        else:
                            sf_phase = cachedSolTabs[solTabList[0]._v_title]
                        sf_phase.setSelection(ant=[ant], pol=[pol], dir=[dir], freq=freqs.tolist(),
                            time={'min':np.min(times), 'max':np.max(times)})
                        val_amp = sf.getValues()[0]
                        val_phase = sf_phase.getValues()[0]
                        val = val_amp * np.cos(val_phase)
                    elif parm == 'Imag':
                        solTabList = getSoltabFromSolType(solType, solTabs, parm='ampl')
                        if not solTabList[0]._v_title in cachedSolTabs:
                            sf_amp = solFetcher(solTabList[0], useCache=True)
                            cachedSolTabs[solTabList[0]._v_title] = sf_amp
                        else:
                            sf_amp = cachedSolTabs[solTabList[0]._v_title]
                        sf_amp.setSelection(ant=[ant], pol=[pol], dir=[dir], freq=freqs.tolist(),
                            time={'min':np.min(times), 'max':np.max(times)})
                        val_phase = sf.getValues()[0]
                        val_amp = sf_amp.getValues()[0]
                        val = val_amp * np.sin(val_phase)
                    else:
                        val = sf.getValues()[0]

                    # Apply flags
                    weights = sf.getValues(weight=True)[0]

                    # etienne part; if it is borken, curse his name
                    # check whether this is clock or tec; if so, reshape properly to account for all freqs in the parmdb
                    # anyway these tables are freq-indep
                    #if solType == "Clock" or solType == "TEC" or solType == "RotationMeasure":
                    #    # find freq-dimensionality
                    #    nfreq = freqs.shape[0]
                    #    print val.shape
                    #    # reshape such that all freq arrays are filled properly
                    #    val = np.tile( val, np.append([nfreq], np.ones(len(val.shape)) ) )
                    #    print val.shape
                    #    weights = np.tile( weights, np.append([nfreq], np.ones(len(weights.shape)) ) )

                    flags = np.zeros(shape=weights.shape, dtype=bool)
                    flags[np.where(weights == 0)] = True
                    if parm == 'Real':
                        weights2 = sf_phase.getValues(weight=True)[0]
                        flags[np.where(weights2 == 0)] = True
                    if parm == 'Imag':
                        weights2 = sf_amp.getValues(weight=True)[0]
                        flags[np.where(weights2 == 0)] = True
                    np.putmask(val, flags, np.nan)

                    shape = data_out[solEntry]['values'].shape
                    #print "shape"
                    #print 'parmdb', shape
                    #print 'h5parm', val.shape
                    #print sf.getAxesNames()
                    #print sf.getType()
                    #print "parmdb", times
                    #for t in times: print '%.1f' % t
                    #print "h5parm", sf.time
                    #for t in sf.time: print '%.1f' % t
                    try:
                        data_out[solEntry]['values'] = val.T.reshape(shape)
                    except ValueError, err:
                        logging.critical('Mismatch between parmdb table and H5parm '
                        'solution table: Differing number of frequencies and/or times')
                        sys.exit(1)
                ipbar += 1
                pbar.update(ipbar)
            else:
                # Handle TECScreen parmdb
                #
                # Get timewidths, freqwidth and freq from first (non-TEC, phase)
                # solentry
                for nonTECsolType in pdbSolTypes:
                    if nonTECsolType != 'TECScreen' and 'Phase' in nonTECsolType:
                        break
                parmname = pdb_in.getNames(nonTECsolType+':*')[0]
                timewidths = pdb_in.getValuesGrid(parmname)[parmname]['timewidths']
                freqwidth = pdb.getValuesGrid(parmname)[parmname]['freqwidths'][0]
                freq = pdb.getValuesGrid(parmname)[parmname]['freqs'][0]
                data_out = makeTECparmdb(h5parm_in, solset, st_tec, timewidths, freq, freqwidth)

            pdb_out.addValues(data_out)

        pbar.finish()

    logging.info('Done.')
