"""
multirefpredict.py
Automated workflow to predict multireference character of molecules in quantum chemistry calculation

Handles the primary functions
"""
from abc import ABC, abstractmethod

class Diagnostic(ABC):
    @abstractmethod
    def computeDiagnostic(self):
        pass
