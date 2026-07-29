import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from time import perf_counter
from src.system import System

def compare_force_methods(n_cell_list, n_repeats):
    results=[]
    for n_cell in n_cell_list:
        s=System(n_cell, rho_star=0.844, dt=0.0001, r_c=2.5, T_star=0.71)
        times_v1=[]
        for repeat in range(n_repeats):
            start=perf_counter()
            s.compute_forces_v1()
            end=perf_counter()
            times_v1.append(end-start)
        times_v1=min(times_v1)
        times_v2=[]
        for repeat in range(n_repeats):
            start=perf_counter()
            s.compute_forces_v2()
            end=perf_counter()
            times_v2.append(end-start)
        times_v2=min(times_v2)
        difference = times_v1 - times_v2
        results.append((n_cell, s.N, times_v1, times_v2, difference))
    return results

results = compare_force_methods([2,3,4,5,6,7], 5)
for r in results:
    print(r)