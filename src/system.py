import numpy as np
class System:
    def __init__(self, n_cell, rho_star, dt, r_c, T_star): #creates FCC lattice at target density
        self.n_cell = n_cell
        self.rho_star = rho_star
        self.dt = dt
        self.r_c = r_c
        self.T_star = T_star
        self.N = 4 * n_cell**3
        self.L_star = (self.N / rho_star)**(1/3)
        offsets=[(0, 0, 0), (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5)]
        a=self.L_star/self.n_cell
        positions_list=[]
        for i in range(self.n_cell):
            for j in range(self.n_cell):
                for k in range(self.n_cell):
                    corner = (i*a, j*a, k*a)
                    for offset in offsets:
                        x = corner[0] + offset[0]*a
                        y = corner[1] + offset[1]*a
                        z = corner[2] + offset[2]*a
                        positions_list.append((x, y, z))
        self.positions = np.array(positions_list)
        velocities=np.random.normal(0, np.sqrt(T_star), (self.N, 3))
        velocities -= np.mean(velocities, axis=0)
        self.velocities = velocities
    def compute_forces(self):
        forces = np.zeros((self.N, 3))
        potential_energy = 0
        U_shift = 4*((1/self.r_c)**12 - (1/self.r_c)**6)
        for i in range(self.N):
            for j in range(i+1, self.N):
                displacement_ij = self.positions[i] - self.positions[j]
                displacement_ij -= np.round(displacement_ij / self.L_star) * self.L_star
                distance_ij=np.linalg.norm(displacement_ij)
                if distance_ij<self.r_c:
                    F_scalar=24/distance_ij*(2*(1/distance_ij)**12-(1/distance_ij)**6)
                    direction = displacement_ij / distance_ij
                    F_vector = F_scalar * direction
                    forces[i] += F_vector
                    forces[j] -= F_vector
                    potential_energy += 4*((1/distance_ij)**12 - (1/distance_ij)**6) - U_shift
        self.forces=forces
        self.potential_energy=potential_energy
        return self.forces, self.potential_energy
    def step(self):
        accelerations = self.forces.copy()
        self.positions += self.velocities * self.dt + 0.5 * accelerations * self.dt**2
        self.compute_forces()
        avg_accelerations = (self.forces + accelerations) / 2
        self.velocities += avg_accelerations * self.dt
        self.positions %= self.L_star
    def compute_temperature(self):
        kinetic_energy = 0.5 * np.sum(self.velocities**2)
        T_inst = (2 * kinetic_energy) / (3 * self.N)
        self.T_inst = T_inst
        self.kinetic_energy = kinetic_energy
        return self.T_inst, self.kinetic_energy
    def run(self, n_production_steps, sample_interval):
        # equilibration variables
        total_energies=[]
        temperatures=[]
        step_count=0
        equilibration=False
        # production variables
        n_snapshots = n_production_steps // sample_interval
        trajectories = np.zeros((n_snapshots, self.N, 3))
        prod_temps=[]
        prod_Enrgs=[]
        while (equilibration == False):
            self.step()
            self.compute_forces()
            self.compute_temperature()
            self.velocities*=np.sqrt(self.T_star/self.T_inst) # isokinetic scaling
            total_energies.append(self.kinetic_energy + self.potential_energy)
            temperatures.append(self.T_inst)
            step_count+=1
            if (step_count % 1000 == 0 and step_count>=2000):
                current_1000_temps = temperatures[-1000:]
                current_temp_1000_avg = np.mean(current_1000_temps)
                current_1000_enrg = total_energies[-1000:]
                previous_1000_enrg = total_energies[-2000:-1000]
                current_enrg_1000_avg = np.mean(current_1000_enrg)
                previous_enrg_1000_avg = np.mean(previous_1000_enrg)
                if (abs(current_temp_1000_avg - self.T_star)/self.T_star<0.005 and abs(current_enrg_1000_avg - previous_enrg_1000_avg)/previous_enrg_1000_avg<0.005):
                    equilibration=True
                    print(f"Equilibration achieved at step {step_count}.")
            if (step_count > 20000):
                print("Equilibration not achieved within 20000 steps.")
                break
        for i in range(n_production_steps):
            self.step()
            self.compute_forces()
            self.compute_temperature()
            if (i % sample_interval == 0):
                slot = i // sample_interval
                trajectories[slot] = self.positions.copy()
                prod_temps.append(self.T_inst)
                prod_Enrgs.append(self.kinetic_energy+self.potential_energy)
        return trajectories, prod_temps, prod_Enrgs