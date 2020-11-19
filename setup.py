# -*- coding: utf-8 -*-

# Learn more: https://github.com/kennethreitz/setup.py

from setuptools import setup, find_packages


with open('README.rst') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='opposed_premix',
    version='0.1.0',
    description='A package for Combustion simulation in a geometory of 1D opposed premix flow.',
    long_description=readme,
    author='Akira Shioyoke',
    author_email='s.akira2986@gmail.com',
    install_requires=['cantera', 'numpy', 'pandas', 'tqdm', 'scipy', 'matplotlib'],
    url='',
    license=license,
    packages=find_packages(exclude=('tests', 'docs'))
)

