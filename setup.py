# Use pip install -e .

from setuptools import setup, find_packages

setup(
    name='CAMPAREE',
    version='0.3',
    packages=find_packages(),
    scripts=["bin/run_camparee.py"],
)
