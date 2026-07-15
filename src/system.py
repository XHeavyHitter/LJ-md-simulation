import numpy as np
class System:
    def __init__(self, n_cell, rho_star, dt, r_c, T_star):
        #atributes
        self.n_cell = n_cell
        self.rho_star = rho_star
        self.dt = dt
        self.r_c = r_c
        self.T_star = T_star
        self.N = 4 * n_cell**3
        #positions
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
        #velocities
        velocities=np.random.normal(0, np.sqrt(T_star), (self.N, 3))
        velocities -= np.mean(velocities, axis=0)
        self.velocities = velocities
    
    def compute_forces(self):
        forces = np.zeros((self.N, 3))
        potential_energy = 0
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
                    potential_energy += 4*((1/distance_ij)**12 - (1/distance_ij)**6)
        self.forces=forces
        self.potential_energy=potential_energy
        return self.forces, self.potential_energy
         