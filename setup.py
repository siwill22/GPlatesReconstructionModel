#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages
#import versioneer

#versioneer.versionfile_source = 'GPlatesReconstructionModel/_version.py'
#versioneer.versionfile_build = 'GPlatesReconstructionModel/_version.py'
#versioneer.tag_prefix = ''
#versioneer.parentdir_prefix = 'GPlatesReconstructionModel-'


install_requires = ['numpy', 'scipy', 'pandas', 'geopandas', 'matplotlib', 'xarray', 'pygmt', 'pooch', 'rasterio', 'PlateTectonicTools']

setup(name='GPlatesReconstructionModel',
      #version=versioneer.get_version(),
      #cmdclass=versioneer.get_cmdclass(),
      description='GPlates ',
      #long_description=long_description,
      url='https://github.com/siwill22/GPlatesReconstructionModel',
      authors='Simon Williams',
      #author_email='siwill22@gmail.com',
      #license='MIT',
      classifiers=[
          'Intended Audience :: Science/Research',
          'Intended Audience :: Developers',
          #'License :: OSI Approved :: BSD License',
          'Natural Language :: English',
          'Operating System :: OS Independent',
          'Programming Language :: Python',
          'Programming Language :: Python :: 3',
          'Topic :: Scientific/Engineering'
      ],
      keywords=[''],
      packages=find_packages(),
      include_package_data=True,
      install_requires=install_requires,
      python_requires='>=3.8')