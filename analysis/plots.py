import matplotlib.pyplot as plt
import numpy as np
def energy_vs_time (total_energies, prod_Enrg, step_count, sample_interval, dt):
    equilibration_times = np.arange(len(total_energies)) * dt
    production_times = step_count * dt + np.arange(len(prod_Enrg)) * sample_interval * dt
    time = np.concatenate([equilibration_times, production_times])
    Energies = total_energies + prod_Enrg
    plt.plot(time, Energies)
    plt.title("Total energy vs time")
    plt.xlabel("Time")
    plt.ylabel("Total energy")
    equilibration_t = step_count * dt
    plt.vlines(equilibration_t, ymin=min(Energies), ymax=max(Energies), colors='red', linestyles='dashed', label='Equilibration end')
    plt.legend()
    plt.show()
def rdf_plot (r_values, mean_g, std_g):
    plt.errorbar(r_values, mean_g, yerr=std_g, fmt='-')
    plt.title("Radial distribution function")
    plt.xlabel("r")
    plt.ylabel("g(r)")
    plt.hlines(1, xmin=min(r_values), xmax=max(r_values), colors='red', linestyles='dashed', label='Reference axis')
    plt.legend()
    plt.show()