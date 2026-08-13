import numpy as np
def rdf(trajectories, L_star, n_blocks=8, bin_width=1/50, max_radius=None):
    if max_radius is None:
        max_radius = L_star / 2
    trajectories = trajectories % L_star # PBC applied
    N=trajectories.shape[1]
    rho_star=N/L_star**3
    bins=int(max_radius/bin_width)
    bin_edges=np.linspace(0, max_radius, bins+1)
    blocks = np.array_split(trajectories, n_blocks)
    r_low = bin_edges[:-1]   # all edges except the last — the "low" side of each bin
    r_high = bin_edges[1:]   # all edges except the first — the "high" side of each bin
    V_shell = (4/3) * np.pi * (r_high**3 - r_low**3)
    all_block_g = [] # holds the g(r) for each block
    for block in blocks:
        block_histogram = np.zeros(bins) # how many particle pairs fall in each bin, summed across all snapshots in this block
        for snapshot in block:
            diff=snapshot[np.newaxis, :, :]-snapshot[:, np.newaxis, :]
            diff=diff-L_star*np.round(diff/L_star) #Applying the minimum image convention
            distances=np.linalg.norm(diff, axis=2) #Scalar distance for every i,j pair
            i_upper, j_upper = np.triu_indices(N, k=1)  # k=1 skips the diagonal
            unique_distances = distances[i_upper, j_upper]  # 1D array of unique pairwise distances
            counts, _ = np.histogram(unique_distances, bins=bin_edges)
            counts/=len(block) # average over snapshots in the block
            block_histogram += counts
        normalised_g=block_histogram/(N*rho_star*V_shell/2) #divided by expected count in a uniform gas
        all_block_g.append(normalised_g)
    all_block_g = np.array(all_block_g)  # shape (n_blocks, bins)
    mean_g = np.mean(all_block_g, axis=0)
    std_g = np.std(all_block_g, axis=0)
    r_values = (bin_edges[:-1] + bin_edges[1:]) / 2
    return r_values, mean_g, std_g