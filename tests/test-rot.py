from ..UB import rot
import numpy as np
r30x = np.array([[1., 0., 0.],
                 [0., 0.8660254, -0.5],
                 [0., 0.5, 0.8660254]])
r30y = np.array([[0.8660254, 0, -0.5],
                 [0., 1., 0.],
                 [0.5, 0., 0.8660254]])
r30z = np.array([[0.8660254, -0.5, 0.],
                [0.5, 0.8660254, 0.],
                [0., 0., 1.]])

assert ((rot(30, 'x') - r30x) < 1e-5).all()
assert ((rot(30, 'y') - r30y) < 1e-5).all()
assert ((rot(30, 'z') - r30z) < 1e-5).all()
