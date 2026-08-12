import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'src'))
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'analysis'))
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'benchmarks'))
from force_methods import compare_force_methods
from plots import benchmark_plot, benchmark_table
results = compare_force_methods(n_cell_list=[2, 3, 4, 5, 6, 7], n_repeats=5)
benchmark_plot(results)
benchmark_table(results)