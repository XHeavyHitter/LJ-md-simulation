import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '../src'))
from system import System
import numpy as np
# Force calculation method version is not fixed, it can be changed to v1, v2, or v3 depending on which method you want to test.
# In all cases we consider a 2 particle system to test the methods against hand calculated values
def test_force_symmetry(): # Checks Newton's 3rd law
    s=System(n_cell=1, rho_star=0.844, dt=0.001, r_c=2.5, T_star=0.71)
    s.positions=np.array([[0, 0, 0], [1.5, 0, 0]])
    s.N=2
    s.compute_forces_v3()
    assert np.allclose(s.forces[0], -s.forces[1])
def test_energy_calculation(): # Checks that the potential energy is calculated correctly
    s=System(n_cell=1, rho_star=0.844, dt=0.001, r_c=2.5, T_star=0.71)
    s.positions=np.array([[0, 0, 0], [1.5, 0, 0]])
    s.N=2
    s.L_star = 10 # Default L_star at n_cell=1 is small enough to cause wrapping
    U_shift = 4*((1/s.r_c)**12 - (1/s.r_c)**6)
    s.compute_forces_v3()
    expected_potential_energy = 4*((1/1.5)**12 - (1/1.5)**6) - U_shift
    assert np.isclose(s.potential_energy, expected_potential_energy)
def test_PBC(): # Checks that the minimum image convention is applied correctly
    s=System(n_cell=1, rho_star=0.844, dt=0.001, r_c=2.5, T_star=0.71)
    s.N=2
    s.positions=np.array([[0, 0, 0], [s.L_star-0.5, 0, 0]]) # True minimum image distance is 0.5
    s.compute_forces_v3()
    F_scalar=24/0.5*(2*(1/0.5)**12-(1/0.5)**6)
    assert np.isclose(np.linalg.norm(s.forces[0]), F_scalar)
