
from setuptools import setup, find_packages
from os import path

# Package meta-data.
NAME = 'mppi'
DESCRIPTION = 'Multi Purpose Python Interface. A python package for managing computations and performing \
post-processing in QuantumESPRESSO and Yambo'
URL = 'https://github.com/marcodalessandro76/MPPI'
EMAIL = 'marco.dalessandro@ism.cnr.it'
AUTHOR = "Marco D'Alessandro"
REQUIRES_PYTHON = '>=3.8' # f-strings and subprocess.run(capture_output=...) are used
VERSION = '1.3'
REQUIRED = ['numpy>=1.20','scipy>=1.6','matplotlib>=3.3','netCDF4>=1.5']
EXTRAS = {
    'docs' : ['sphinx','sphinx_rtd_theme','nbsphinx','ipykernel'], # to build the documentation
    'test' : ['pytest']                                             # to run the test suite
}
#-----------------------------------------------

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name=NAME,
    version=VERSION,
    description=DESCRIPTION,
    long_description=long_description,
    long_description_content_type='text/markdown',
    url=URL,
    author=AUTHOR,
    author_email=EMAIL,
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Physics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3 :: Only'
    ],
    keywords='python system-interface post-processing QuantumESPRESSO Yambo',
    packages=find_packages(exclude=['contrib', 'docs', 'tests']),
    python_requires=REQUIRES_PYTHON,
    install_requires=REQUIRED,
    extras_require=EXTRAS
)
