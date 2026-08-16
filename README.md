# LJ-md-simulation
Molecular dynamics simulation using the Lennard-Jones potential, with RDF validation against literature data
## Physics background
<p>Lennard-Jones potential energy:<br>
$$U(r)=4\epsilon \left[\left(\frac{\sigma}{r}\right)^{12}-\left(\frac{\sigma}{r}\right)^6\right]$$

<p>LJ potential models atomic interactions: $r^{-12}$ term (Pauli repulsion) + $r^{-6}$ term (van der Waals attraction). Simulation runs in reduced units ($\sigma=\epsilon=1$)</p> 

<p>Force equation:<br>

<p>$$F(r)=\frac{24\epsilon}{r}\left[2\left(\frac{\sigma}{r}\right)^{12}-\left(\frac{\sigma}{r}\right)^6\right]$$<br>

Full derivation: [fundamentals.2_particle_test.ipynb](notebooks/fundamentals.2_particle_test.ipynb)
## Results
<p>Radial distribution function of a Lennard-Jones fluid was validated near the triple point ($\rho*=0.844$, $T*=0.71$, $N=500$).<br> 
![rdf_plot.png](results/rdf_plot.png)

The obtained plot qualitatively matches the ones in Frenkel and Smit's textbook, and Rahman's paper. The first RDF peak has a value $g(r^*)=3.00\pm 0.06$, $r^*=1.10$. Because Rahman's paper and Frenkel and Smit's textbook compute their respective plots using slightly different parameters, the comparison between obtained values and literature values, in this case, is unfair. Nevertheles, since the obtained values are in the ballpark of the literature ones and as all of the parameters are considered to be in the triple point region, it can be concluded that the obtained $g(r^*)$ is correct. The obtained diffusion constant ($D^*$) is also within the comparable order of magnitude of the one that is obtained in Rahman's paper. After numerous runs the obtained value is $D^*=0.021\pm 0.004$. 
## Method
| Parameter | Value |
|---|---|
| N | $500$ |
| $T^*$ | $0.71$ |
| $\rho^*$ | $0.844$ |
| $r_c$ | $2.5$ |
| $\Delta t^*$ | $0.0001$ |
| $L^*$ | $8.40$ |
| Steps | $50000$ |

<p>The force calculation method has three versions.<br> 
![force_methods_benchmark.png](results/force_methods_benchmark.png) ![force_methods_benchmark_table.png](results/force_methods_benchmark_table.png)
 
The results suggest that version 3 is the fastest at 500 particles, with version 2 being a close second, and a faster method at smaller quantities of particles.
## References
1. Rahman, A. Correlations in the Motion of Atoms in Liquid Argon. Phys. Rev. 1964, 136 (2A), A405–A411.
2. Frenkel, D.; Smit, B. Understanding Molecular Simulation: From Algorithms to Applications, 2nd ed.; Academic Press: San Diego, 2002.