import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '../src'))
import pytest
import numpy as np
from system import System

@pytest.fixture(scope="module")
def conservation_data():
    s = System(n_cell=2, rho_star=0.844, dt=0.0001, r_c=2.5, T_star=0.71)
    s.compute_forces_v3()
    total_energy = []
    momentum = []
    for i in range(1000):
        s.step()
        s.compute_temperature()
        total_energy.append(s.potential_energy + s.kinetic_energy)
        momentum.append(np.sum(s.velocities, axis=0))
    return total_energy, momentum