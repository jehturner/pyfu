""" Add the pyfu directory to start of the Python path
"""
import iraf,sys

pyfu_path = iraf.osfn('pyfu$')

if pyfu_path not in sys.path:
    sys.path.insert(1,pyfu_path)

