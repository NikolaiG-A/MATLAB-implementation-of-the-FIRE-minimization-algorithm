# MATLAB-implementation-of-the-FIRE-minimization-algorithm

These files contain MATLAB implementation of the fast inertial relaxation engine (FIRE) algorithm. It is used in local energy minimisation during atomistic simulations in Molecular Dynamics (MD) simulations. The implementation follows the original article Guénolé, J., Nöhring, W. G., Vaid, A., Houllé, F., Xie, Z., Prakash, A., & Bitzek, E. (2020). Assessment and optimization of the fast inertial relaxation engine (fire) for energy minimization in atomistic simulations and its implementation in lammps. Computational Materials Science, 175, 109584.

The files include:
1. Parameter initialization (params.m).
2. The MD integration schemes including the explicit Euler integration, semi-implicite Euler integration, Leapfrog integration, Verlet integration (MDintegrator.m). The integration algorithms are adjusted as described in the paper.
3. The FIRE minimization algorithm (FIRE.m).
4. Example file (example.m) of minimization of the multivalued function
<img src="https://render.githubusercontent.com/render/math?math= f = e^{x-2*x^2-y^2}sin(6*(x + y + x*y^2))">
