LoSoTo: LOFAR solutions tool

Author:
Francesco de Gasperin (overall project)
Cyril Tasse (APPLY)
Reinout Van Weeren (SMOOTH, CLOCKTEC)
Maaijke Mevius (CLOCKTEC)
Bas van der Tol (CLOCKTEC, TOOLS)

This library have tools to deal with LOFAR solution tables.

- set to zero solutions (step: RESET)
- clock/tec separation (step: CLOCKTEC)
- flag outliers (step: FLAG)
- solutions smoothing (step: SMOOTH)
- solutions interpolation (step: INTERP)
- make some inspection plots (step: PLOT)
- apply solutions (step: APPLY)

Included tools:

- H5parm_creator.py:
    * creates an h5parm file from an MS or a globaldb created with parmdb_collector.py
- H5parm_merge.py:
    * merge multiple H5parm files into a single one
- parmdb_collector.py:
    * fetches parmdb tables from the cluster

Required packages:

- pytables 3.0.0 (numpy >= 1.4.1 & numexpr >= 2.0.0)
