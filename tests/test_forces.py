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
def test_energy_calculation():
    s=System(n_cell=1, rho_star=0.844, dt=0.001, r_c=2.5, T_star=0.71)
    s.positions=np.array([[0, 0, 0], [1.5, 0, 0]])
    s.N=2
    s.L_star = 10
    U_shift = 4*((1/s.r_c)**12 - (1/s.r_c)**6)
    s.compute_forces()
    expected_potential_energy = 4*((1/1.5)**12 - (1/1.5)**6) - U_shift
    assert np.isclose(s.potential_energy, expected_potential_energy)
def test_PBC():
    s=System(n_cell=1, rho_star=0.844, dt=0.001, r_c=2.5, T_star=0.71)
    s.N=2
    s.positions=np.array([[0, 0, 0], [s.L_star-0.5, 0, 0]])
    s.compute_forces()
    F_scalar=24/0.5*(2*(1/0.5)**12-(1/0.5)**6)
    assert np.isclose(np.linalg.norm(s.forces[0]), F_scalar)
