import os
import re
import io

from setuptools import find_packages, setup

with open(os.path.join(os.path.dirname(__file__), 'README.md')) as readme:
    README = readme.read()
# allow setup.py to be run from any path
os.chdir(os.path.normpath(os.path.join(os.path.abspath(__file__), os.pardir)))


setup(
    name='flutype_analysis',
    packages=find_packages(),
    include_package_data=True,
    license='LGPLv3',
    description='Analysis module for virus-peptide interactions',
    long_description=README,
    author='Jan Grzegorzewski & Matthias König',
    author_email='janekg89@hotmail.de',
)