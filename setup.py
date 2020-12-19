#!/usr/bin/env python
#-----------------------------------------------------------------------------
# Copyright (c) 2020, Greg Landrum
# All rights reserved.
#
# The full license is in the LICENSE file, distributed with this software.
#-----------------------------------------------------------------------------

from setuptools import setup, find_packages

INSTALL_REQUIRES = ['intake >=0.5.2']

setup(
    name='intake-rdkit',
    version="0.1.0",
    description='rdkit plugins for Intake',
    url='https://github.com/greglandrum/intake-rdkit',
    maintainer='greg landrum',
    maintainer_email='greg.landrum@t5informatics.com',
    license='BSD',
    py_modules=['intake_rdkit'],
    packages=find_packages(),
    entry_points={
        'intake.drivers': [
            'sdf = intake_rdkit.sdf:SDFSource',
        ]
    },
    package_data={'': ['*.csv', '*.yml', '*.html']},
    include_package_data=True,
    install_requires=INSTALL_REQUIRES,
    long_description=open('README.md').read(),
    long_description_content_type="text/markdown",
    zip_safe=False, )
