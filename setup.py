#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import warnings

from setuptools import setup, Extension, find_packages

install_requires = [
                    'tomopy>=0.0.3'
                    ]

# Main setup configuration.
setup(
      name='mapstotomopy',
      version=open('VERSION').read().strip(),
      
      packages = ['mapstotomopy'],
      include_package_data = True,
      
      author='Young Hong',
      author_email='younghong2015@u.northwestern.edu',

      description='use h5 files from maps as an input',
      keywords=['tomography', 'reconstruction', 'imaging'],
 

      license='BSD',
      platforms='Any',
      install_requires = install_requires,

      classifiers=['Development Status :: 4 - Beta',
		   'License :: OSI Approved :: BSD License',
		   'Intended Audience :: Science/Research',
		   'Intended Audience :: Education',
		   'Intended Audience :: Developers',
		   'Natural Language :: English',
		   'Operating System :: OS Independent',
		   'Programming Language :: Python',
		   'Programming Language :: Python :: 2.6',
		   'Programming Language :: Python :: 2.7',
		   'Programming Language :: C',
		   'Programming Language :: C++']
      )
