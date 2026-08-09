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
        times_v3=[]
        for repeat in range(n_repeats):
            start=perf_counter()
            s.compute_forces_v3()
            end=perf_counter()
            times_v3.append(end-start)
        times_v3=min(times_v3)
        difference12=abs(times_v1-times_v2)
        difference13=abs(times_v1-times_v3)
        difference23=abs(times_v2-times_v3)
        results.append((n_cell, s.N, times_v1, times_v2, times_v3, difference12, difference13, difference23))
    return results