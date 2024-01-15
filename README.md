# Monte Carlo Method
## Reference
A. L. Garcia, "Numerical Methods for Physics": \
 Chapter 11 Stochastic Methods
## Exercisies
### 11-2.dsmceq
Modify  **dsmceq.cpp**  to compute the expected equilibrium speed distribution histogram using the Maxwell-Boltzmann distribution (see Exercise 11.7). Plot this distribution along with the measured speed distribution histogram and demonstrate that the program correctly approaches equilibrium.

### 11-3.dsmcne_color
Modify  **dsmcne.cpp**  so that particles are labeled “black” or “white”. Particles that reflect off a wall are turned black with probability $q$ and white with probability $1-q$. Set the walls at equal temperature and make them stationary; give them different values of $q$ to set up a pigmentation gradient across the system (see Figure 11.10). Measure the average pigment flux, and compute the self-diffusion coefficient using (11.54). Compare your results with
```math
D=\frac{6\pi}{32} \left \langle v \right \rangle \lambda
```
the value given by Chapman-Enskog theory.
