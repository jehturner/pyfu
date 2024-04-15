""" Add the pyfu directory to start of the Python path
"""
import sys, os.path
import iraf

pyfu_dir = iraf.osfn('pyfu$')

parent_dir, pkg = os.path.split(os.path.normpath(pyfu_dir))
if pkg != 'pyfu':
    raise ValueError('Path pyfu$ in IRAF must end in "/pyfu/')

if parent_dir not in sys.path:
    sys.path.insert(1, parent_dir)

