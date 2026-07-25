import matplotlib.pyplot as plt
import numpy as np
def energy_vs_time(total_energies, prod_Enrg, step_count, sample_interval, dt):
    equilibration_times = np.arange(len(total_energies)) * dt
    production_times = step_count * dt + np.arange(len(prod_Enrg)) * sample_interval * dt
    time = np.concatenate([equilibration_times, production_times])
    Energies = total_energies + prod_Enrg
    plt.plot(time, Energies)
    plt.title("Total energy vs time")
    plt.xlabel("Time")
    plt.ylabel("Total energy")
    equilibration_t = step_count * dt
    plt.vlines(equilibration_t, ymin=min(Energies), ymax=max(Energies), colors='red', linestyles='dashed')
    plt.show()