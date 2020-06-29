from ..UB import Bmat_gen
import numpy as np

Bmatck = np.eye(3)
Bmat = Bmat_gen(np.array([2*np.pi, 2*np.pi, 2*np.pi]), np.array([90, 90, 90]))
assert (np.abs(Bmat - Bmatck) < 1e-5).all()
