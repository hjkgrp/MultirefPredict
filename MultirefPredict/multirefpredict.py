"""
multirefpredict.py
Automated workflow to predict multireference character of molecules in quantum chemistry calculation

Handles the primary functions
"""
from MultirefPredict.b1 import B1

def diagnostic_factory(diagnostic_type, **kwargs):
    if diagnostic_type == "B1":
        return B1(**kwargs)
    else:
        raise Exception("Diagnostic type not found")
