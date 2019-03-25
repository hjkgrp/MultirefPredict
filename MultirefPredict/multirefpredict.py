"""
multirefpredict.py
Automated workflow to predict multireference character of molecules in quantum chemistry calculation

Handles the primary functions
"""
from MultirefPredict.ebased_diagnostic import B1,A25PBE,TAE
from MultirefPredict.ccbased_diagnostic import CCBased
from MultirefPredict.casbased_diagnostic import C0

def diagnostic_factory(diagnostic_type, **kwargs):
    """
    factory class for calculating diagnostics

    Parameters
    ----------
    diagnostic_type : str
        currently support: {"B1","A25PBE"}
        The type of diagnostic to calculate
    keyword arguments:
        molecule: qcelement molecule instance

    Returns
    -------
    cls_instance:  instance of a class for calculating specific diagnostics
    """
    cls_dict = dict(B1=B1,A25PBE=A25PBE,TAE=TAE, 
                    T1=CCBased, D1=CCBased, D2=CCBased,
                    CCBased=CCBased,
                    C0=C0)

    if diagnostic_type not in cls_dict.keys():
        raise Exception("Diagnostic type not found")

    cls = cls_dict[diagnostic_type]
    cls_instance = cls(**kwargs)
    return cls_instance
