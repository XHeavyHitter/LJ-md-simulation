import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'src'))
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'analysis'))
from system import System
from plots import energy_vs_time, energy_vs_time_production, temperature_vs_time, rdf_plot
from radial_distribution_function import rdf
from mean_square_displacement import msd, diffusion
system=System(n_cell=5, rho_star=0.844, T_star=0.71, dt=0.0001, r_c=2.5)
trajectories, prod_temps, prod_Enrgs, total_energies, temperatures, step_count = system.run(n_production_steps=50000, sample_interval=50)
r_values, mean_g, std_g, r_peak1, g_peak1, g_peak1_std = rdf(trajectories, system.L_star)
MSD=msd(trajectories)
D=diffusion(MSD, system.sample_interval, system.dt)
energy_vs_time(total_energies, prod_Enrgs, step_count, system.sample_interval, system.dt)
energy_vs_time_production(prod_Enrgs, step_count, system.sample_interval, system.dt)
temperature_vs_time(temperatures, prod_temps, step_count, system.sample_interval, system.dt)
rdf_plot(r_values, mean_g, std_g)
with open('results/simulation_summary.txt', 'w') as f:
    f.write(f"Diffusion coefficient D: {D}\n")
    f.write("\n")
    f.write(f"RDF first peak position r*: {r_peak1}\n")
    f.write(f"RDF first peak height g(r*): {g_peak1} +/- {g_peak1_std}\n")