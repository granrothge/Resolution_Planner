import numpy as np


class Disc_chop(object):
    """
    a class to define a disc chopper for resolution calculation
    """

    def __init__(self, R, dphi):
        self.R = R
        self.dphi = dphi

    def Deltat(self, f=60.0):
        """
        f is the frequency in Hz
        returns the dt
        """
        return self.dphi/(2*np.pi*f)
