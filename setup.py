"""
MultirefPredict
Automated workflow to predict multireference character of molecules in quantum chemistry calculation
"""
from setuptools import setup, find_packages
import versioneer

short_description = __doc__.split("\n")

try:
    with open("README.md", "r") as handle:
        long_description = handle.read()
except:
    long_description = "\n".join(short_description[2:]),


setup(
    # Self-descriptive entries which should always be present
    name='MultirefPredict',
    author='Kulik group at MIT',
    author_email='fangliu@mit.edu',
    description=short_description[0],
    long_description=long_description,
    long_description_content_type="text/markdown",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    license='BSD-3-Clause',

    # Which Python importable modules should be included when your package is installed
    #packages=['MultirefPredict', "MultirefPredict.tests"],
    packages=find_packages(),

    # Optional include package data to ship with your package
    # Comment out this line to prevent the files from being packaged with your software
    # Extend/modify the list to include/exclude other items as need be
    package_data={'MultirefPredict': ["data/*.dat"]
                  },

    # Additional entries you may want simply uncomment the lines you want and fill in the data
    # author_email='me@place.org',      # Author email
    # url='http://www.my_package.com',  # Website
    install_requires=[
        'qcelemental>=0.2.6',
        'qcenegine',
        'openbabel>=2.4.1'
    ],              # Required packages, pulls from pip if needed; do not use for Conda deployment
    tests_require=[
        'pytest',
        'pytest-cov',
    ], 
    # platforms=['Linux',
    #            'Mac OS-X',
    #            'Unix',
    #            'Windows'],            # Valid platforms your code works on, adjust to your flavor
    python_requires=">=3.6",          # Python version restrictions

    # Manual control if final package is compressible or not, set False to prevent the .egg from being made
    # zip_safe=False,

)
