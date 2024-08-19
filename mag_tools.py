import numpy as np


def form_fac(s, A, B, C, D, a, b, c, l=0):
    """ a form factor base on Brown's formaliam
        l is by defulat 0 for higher order terms set this to the approrpiate 
        value.
    """
    jo = A*np.exp(-a*s**2)+B*np.exp(-b*s**2)+C*np.exp(-c*s**2)+D
    if l > 0:
        jo *= s**2
    return jo
