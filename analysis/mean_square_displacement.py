import numpy as np
def msd(trajectories):
    displacement=trajectories-trajectories[0] # per particle, per snapshot displacement from initial lattice
    squared_displacement = np.sum(displacement**2, axis=2)
    msd=np.mean(squared_displacement, axis=1)
    return msd
def diffusion(msd, sample_interval, dt):
    k = np.arange(len(msd))
    time = k * sample_interval * dt
    slope, intercept = np.polyfit(time, msd, 1)
    D = slope / 6
    return D