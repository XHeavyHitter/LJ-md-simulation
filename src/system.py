import numpy as np
import itertools
class System:
    def __init__(self, n_cell, rho_star, dt, r_c, T_star): #creates FCC lattice at target density
        self.n_cell = n_cell
        self.rho_star = rho_star
        self.dt = dt
        self.r_c = r_c
        self.T_star = T_star
        self.N = 4 * n_cell**3
        self.L_star = (self.N / rho_star)**(1/3)
        offsets=[(0, 0, 0), (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5)] # In FCC lattice, each unit cell has 4 atoms at these fractional coordinates
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
        velocities=np.random.normal(0, np.sqrt(T_star), (self.N, 3)) # Velocity components are drawn from a normal distribution with mean 0 and std sqrt(T_star), net following the Maxwell-Boltzmann distribution.
        velocities -= np.mean(velocities, axis=0)
        self.velocities = velocities
    def compute_forces_v1(self): # Base method
        forces = np.zeros((self.N, 3))
        potential_energy = 0
        U_shift = 4*((1/self.r_c)**12 - (1/self.r_c)**6)
        for i in range(self.N):
            for j in range(i+1, self.N):
                displacement_ij = self.positions[i] - self.positions[j] # displacement vector between two particles in the same frame
                displacement_ij -= np.round(displacement_ij / self.L_star) * self.L_star #applying PBC
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
    def compute_forces_v2(self): # Vectorized
         # Force calculation logic
         displacement_ij=self.positions[:, np.newaxis, :]-self.positions[np.newaxis, :, :]
         displacement_ij -= np.round(displacement_ij / self.L_star) * self.L_star #applying PBC
         distance_ij=np.linalg.norm(displacement_ij, axis=2)
         np.fill_diagonal(distance_ij, np.inf) # Makes the diagonal equal to infinities - when values are used, the results are 0
         relevant_pairs=distance_ij<self.r_c
         F_scalar=24/distance_ij*(2*(1/distance_ij)**12-(1/distance_ij)**6)
         F_scalar=F_scalar*relevant_pairs
         direction=displacement_ij/distance_ij[:, :, np.newaxis]
         F_vector = F_scalar[:, :, np.newaxis] * direction
         forces = F_vector.sum(axis=1)
         U_shift = 4*((1/self.r_c)**12 - (1/self.r_c)**6)
         pair_energy = 4*((1/distance_ij)**12 - (1/distance_ij)**6) - U_shift
         pair_energy = pair_energy * relevant_pairs
         i_upper, j_upper = np.triu_indices(self.N, k=1) # The full array double counts every pair and is symmetric
         potential_energy = pair_energy[i_upper, j_upper].sum()
         self.forces = forces
         self.potential_energy = potential_energy
         return self.forces, self.potential_energy
    def compute_forces_v3(self): # Cell lists
         # Cell grid
         n_cells_per_dim = int(self.L_star / self.r_c)
         if n_cells_per_dim <= 0: # Guards against the case where the cutoff radius is larger than the box length.
             n_cells_per_dim = 1
         actual_cell_length = self.L_star / n_cells_per_dim
         cell_indexes = np.array(self.positions / actual_cell_length, dtype=int)
         cell_atoms = {}
         for atom_index in range(self.N):
            cell = tuple(cell_indexes[atom_index])  
            if cell not in cell_atoms:
                cell_atoms[cell] = []
            cell_atoms[cell].append(atom_index)
         # Checking logic (self + 13 forward, avoids double-counting cell-pairs, Newton's third law)
         neighbor_offsets = []
         for offset in itertools.product([-1, 0, 1], repeat=3): # Creates a list of all possible offsets for neighboring cells in 3D
            if offset >= (0, 0, 0):
                neighbor_offsets.append(offset)
         # Force calculation logic
         forces = np.zeros((self.N, 3))
         potential_energy = 0
         for cell in cell_atoms:
            atoms_in_cell = cell_atoms[cell]
            seen_neighbors = set() 
            U_shift = 4*((1/self.r_c)**12 - (1/self.r_c)**6) # Softener for potential energy
            for offset in neighbor_offsets:
                neighbor_cell = tuple((np.array(cell) + np.array(offset)) % n_cells_per_dim)
                if neighbor_cell in seen_neighbors: # Guards against small grids, when multiple cells overlap
                    continue
                seen_neighbors.add(neighbor_cell)
                atoms_in_neighbor = cell_atoms.get(neighbor_cell, [])
                if len(atoms_in_neighbor) == 0:
                    continue  # Nothing in this neighbor cell, skip
                positions_in_cell = self.positions[atoms_in_cell]
                positions_in_neighbor = self.positions[atoms_in_neighbor]
                displacement_ij = positions_in_cell[:, np.newaxis, :] - positions_in_neighbor[np.newaxis, :, :]
                displacement_ij -= np.round(displacement_ij / self.L_star) * self.L_star #applying PBC
                distance_ij = np.linalg.norm(displacement_ij, axis=2)
                if neighbor_cell == cell:
                    np.fill_diagonal(distance_ij, np.inf) # Avoids double coounting for the same cell, samke logic as in v2
                relevant_pairs=distance_ij<self.r_c
                F_scalar=24/distance_ij*(2*(1/distance_ij)**12-(1/distance_ij)**6)
                F_scalar=F_scalar*relevant_pairs
                direction=displacement_ij/distance_ij[:, :, np.newaxis]
                F_vector = F_scalar[:, :, np.newaxis] * direction
                pair_forces = F_vector.sum(axis=1)
                if cell == neighbor_cell:
                    forces[atoms_in_cell] += pair_forces
                else:
                    forces[atoms_in_cell] += pair_forces
                    neighbor_forces = -F_vector.sum(axis=0)
                    forces[atoms_in_neighbor] += neighbor_forces
                # Potential energy calculation
                pair_energy = 4*((1/distance_ij)**12 - (1/distance_ij)**6) - U_shift
                pair_energy = pair_energy * relevant_pairs
                if cell == neighbor_cell:
                    potential_energy += pair_energy.sum() / 2
                else:
                    potential_energy += pair_energy.sum()
         self.forces = forces
         self.potential_energy = potential_energy
         return self.forces, self.potential_energy
    def step(self):
        accelerations = self.forces.copy()
        self.positions += self.velocities * self.dt + 0.5 * accelerations * self.dt**2
        self.compute_forces_v3()
        avg_accelerations = (self.forces + accelerations) / 2
        self.velocities += avg_accelerations * self.dt
    def compute_temperature(self):
        kinetic_energy = 0.5 * np.sum(self.velocities**2)
        T_inst = (2 * kinetic_energy) / (3 * self.N)
        self.T_inst = T_inst
        self.kinetic_energy = kinetic_energy
        return self.T_inst, self.kinetic_energy
    def run(self, n_production_steps, sample_interval):
        self.compute_forces_v3() # Initial force calculation
        # Equilibration variables
        total_energies=[]
        temperatures=[]
        step_count=0
        equilibration=False
        # Production variables
        n_snapshots = n_production_steps // sample_interval
        trajectories = np.zeros((n_snapshots, self.N, 3))
        prod_temps=[]
        prod_Enrgs=[]
        while (equilibration == False): # Equilibration loop; total energy and temperature have to stabilize
            self.step()
            self.compute_temperature()
            self.velocities*=np.sqrt(self.T_star/self.T_inst) # Isokinetic scaling
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
        for i in range(n_production_steps): # Production loop
            self.step()
            self.compute_temperature()
            if (i % sample_interval == 0):
                slot = i // sample_interval
                trajectories[slot] = self.positions.copy()
                prod_temps.append(self.T_inst)
                prod_Enrgs.append(self.kinetic_energy+self.potential_energy)
        return trajectories, prod_temps, prod_Enrgs