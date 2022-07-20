import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

from glob import glob
from os.path import basename
from os.path import splitext

from setuptools import setup
from setuptools import find_packages

def _requires_from_file(filename):
    return open(filename).read().splitlines()

setup(
    name="turbogenius",
    version="0.0.2",
    license="MIT License",
    description="A python-based workflow system for the SISSA quantum Monte Carlo package TurboRVB",
    author="Kousuke Nakano",
    author_email="kousuke_1123@icloud.com",
    url="https://www.kosuke-nakano-research.info",
    entry_points={
        'console_scripts': [
            'turbogenius = turbogenius.turbo_genius_cli:cli',
            'trexio-to-turborvb = turbogenius.trexio_to_turborvb:main',
            ]
    },
    packages=find_packages(),
    include_package_data=True,
    zip_safe=False,
    install_requires=_requires_from_file('requirements.txt')
)
