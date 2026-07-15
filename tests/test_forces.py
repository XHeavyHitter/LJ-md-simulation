import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '../src'))
from system import System
import numpy as np
def test_force_symmetry():
    s=System(n_cell=1, rho_star=0.844, dt=0.001, r_c=2.5, T_star=0.71)
    s.positions=np.array([[0, 0, 0], [1.5, 0, 0]])
    s.N=2
    s.compute_forces()
    assert np.allclose(s.forces[0], -s.forces[1])
