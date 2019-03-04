"""
Helper functions for comparison
"""
import math

def fuzzyEqual(a,b,thre):
    res = False
    if math.fabs(a-b) < thre: 
        res = True
    return res
