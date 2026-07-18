import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '../src'))
from system import System
import numpy as np
def test_temperature():
    s=System(n_cell=1, rho_star=0.844, dt=0.001, r_c=2.5, T_star=0.71)
    s.velocities=np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1], [-1, -1, -1]])
    s.N=4
    s.compute_temperature()
    assert np.isclose(s.T_inst, 0.5)