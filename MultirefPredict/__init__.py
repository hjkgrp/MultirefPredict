"""
MultirefPredict
Automated workflow to predict multireference character of molecules in quantum chemistry calculation
"""

# Add imports here
from .multirefpredict import *

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
