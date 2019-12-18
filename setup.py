#!/usr/bin/env python

from distutils.core import setup

LONG_DESCRIPTION = \
'''
Count statistics of DNA bases in a defined genomic region for over multiple
BAM alignment files.
'''


setup(
    name='base_counter',
    version='0.1.0.0',
    author='Bernie Pope',
    author_email='bjpope@unimelb.edu.au',
    packages=['base_counter'],
    package_dir={'base_counter': 'base_counter'},
    entry_points={
        'console_scripts': ['base_counter = base_counter.base_counter:main',
        'variant_detection = base_counter.variant_detection:main',
        'reference_bases = base_counter.reference_bases:main']
    },
    url='https://github.com/bjpop/base_counter',
    license='LICENSE',
    description=('Count bases in a defined genomic region'),
    long_description=(LONG_DESCRIPTION),
    install_requires=["pysam", "numpy", "scipy", "intervaltree"],
)
