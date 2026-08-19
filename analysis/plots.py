import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
def energy_vs_time (total_energies, prod_Enrg, step_count, sample_interval, dt):
    plt.figure()
    equilibration_times = np.arange(len(total_energies)) * dt
    production_times = step_count * dt + np.arange(len(prod_Enrg)) * sample_interval * dt
    time = np.concatenate([equilibration_times, production_times])
    Energies = total_energies + prod_Enrg
    plt.plot(time, Energies)
    plt.title("Total energy vs time")
    plt.xlabel("t*")
    plt.ylabel("U*")
    equilibration_t = step_count * dt
    plt.vlines(equilibration_t, ymin=min(Energies), ymax=max(Energies), colors='red', linestyles='dashed', label='Equilibration end')
    plt.legend()
    plt.savefig('results/energy_vs_time.png')
def energy_vs_time_production (prod_Enrg, step_count, sample_interval, dt):
    plt.figure()
    t_val=step_count * dt + np.arange(len(prod_Enrg)) * sample_interval * dt
    plt.plot(t_val, prod_Enrg)
    plt.title("Total energy vs time (production phase)")
    plt.xlabel("t*")
    plt.ylabel("U*")
    plt.savefig('results/energy_vs_time_production.png')
def temperature_vs_time (temperatures, prod_temps, step_count, sample_interval, dt):
    plt.figure()
    equilibration_times = np.arange(len(temperatures)) * dt
    production_times = step_count * dt + np.arange(len(prod_temps)) * sample_interval * dt
    time = np.concatenate([equilibration_times, production_times])
    Temps = temperatures + prod_temps
    plt.plot(time, Temps)
    plt.title("Temperature vs time")
    plt.xlabel("t*")
    plt.ylabel("T*")
    equilibration_t = step_count * dt
    plt.vlines(equilibration_t, ymin=min(Temps), ymax=max(Temps), colors='red', linestyles='dashed', label='Equilibration end')
    plt.legend()
    plt.savefig('results/temperature_vs_time.png')
def rdf_plot (r_values, mean_g, std_g):
    plt.figure()
    plt.errorbar(r_values, mean_g, yerr=std_g, fmt='-')
    plt.title("Radial distribution function")
    plt.xlabel("r*")
    plt.ylabel("g(r*)")
    plt.hlines(1, xmin=min(r_values), xmax=max(r_values), colors='red', linestyles='dashed', label='Reference axis')
    plt.legend()
    plt.savefig('results/rdf_plot.png')
def benchmark_plot(results):
    plt.figure()
    N_values = [r[1] for r in results]
    times_v1 = [r[2] for r in results]
    times_v2 = [r[3] for r in results]
    times_v3 = [r[4] for r in results]
    plt.plot(N_values, times_v1, label='compute_forces_v1')
    plt.plot(N_values, times_v2, label='compute_forces_v2')
    plt.plot(N_values, times_v3, label='compute_forces_v3')
    plt.title("Benchmark results")
    plt.xlabel("Number of particles")
    plt.ylabel("Execution time")
    plt.yscale('log')
    plt.legend()
    plt.savefig('results/force_methods_benchmark.png')
def benchmark_table(results):
    plt.figure()
    plt.axis('off')
    table_data = [[r[1], r[5], r[6], r[7]] for r in results]
    column_labels = ["N", "v1-v2 diff", "v1-v3 diff", "v2-v3 diff"]
    plt.table(cellText=table_data, colLabels=column_labels, loc='center')
    plt.savefig('results/force_methods_benchmark_table.png')
def trajectory_animation(trajectories, L_star):
    frames = trajectories % L_star
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.set_xlim(0, L_star)
    ax.set_ylim(0, L_star)
    ax.set_zlim(0, L_star)
    ax.set_box_aspect([1, 1, 1])
    scatter = ax.scatter([], [], [], s=15)
    def update(frame_index):
        positions = frames[frame_index]
        scatter._offsets3d = (positions[:, 0], positions[:, 1], positions[:, 2])
        return scatter,
    ani = FuncAnimation(fig, update, frames=len(frames), interval=50, blit=False)
    ani.save('results/trajectory.mp4', writer='ffmpeg', fps=20)