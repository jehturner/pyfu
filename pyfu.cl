# Copyright(c) 2006 Association of Universities for Research in Astronomy, Inc.

# Add modules to the Python path:
pyexecute("pyfu$addpath.py",verbose=no)

# Define pyifu package tasks:
package pyfu

pyexecute("pyfu$pyfmosaic_iraf.py",tasknames="pyfmosaic,pyfalign")
#pyexecute("pyfu$pyfscatt_iraf.py",tasknames="pyfscatt")

clbye()
