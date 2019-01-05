#!/usr/bin/env python

from setuptools import setup

setup(
    name='cfa_pipeline',
    version='0.0.1',
    author='Bernie Pope, Khalid Mahmood, Jessica Chung, Hongzhi Luo',
    author_email='kmahmood@unimelb.edu.au, hongzhil@unimelb.edu.au',
    packages=['src'],
    entry_points={
        'console_scripts': ['cfapipeline = src.main:main']
    },
    url='https://github.com/verityrise/canfam_pipeline',
    license='LICENSE.txt',
    description='cfapipeline is a pipeline system for bioinformatics workflows\
     with support for running pipeline stages on a distributed compute cluster.',
    long_description=open('README.md').read(),
    install_requires=[
        "ruffus == 2.6.3",
        "drmaa == 0.7.6",
        "pyyaml>=4.2b1"
    ],
)
