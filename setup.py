#!/usr/bin/env python

from distutils.core import setup

setup(name='PyFU',
      version='0.11',
      description='Python scripts for mosaicking/resampling IFU data cubes',
      author='James E.H. Turner',
      author_email='jturner@gemini.edu',
      packages=['pyfu'],
      license='BSD',
      install_requires=[
          'astropy>=3.1',
          'numpy>=1.9,<2',
          'scipy>=0.14',
      ],
      python_requires='>=3.6',
)

