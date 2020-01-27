---
title: 'SARAS: A General-Purpose PDE Solver for Fluid Dynamics'

tags:
  - C++
  - PDE
  - turbulence
  - fluid dynamics

authors:
  - name: Roshan Samuel
    orcid: 0000-0002-1280-9881
    affiliation: 1
  - name: Shashwat Bhattacharya
    orcid: 0000-0001-7462-7680
    affiliation: 1
  - name: Ali Asad
    orcid: 0000-0001-9704-6686
    affiliation: 2
  - name: Soumyadeep Chatterjee
    orcid: 0000-0001-7957-1727
    affiliation: 2
  - name: Mahendra K. Verma
    orcid: 0000-0002-3380-4561
    affiliation: 2

affiliations:
 - name: Department of Mechanical Engineering, Indian Institute of Technology - Kanpur
   index: 1
 - name: Department of Physics, Indian Institute of Technology - Kanpur
   index: 2

date: 15 January 2020

bibliography: paper.bib

---

# Summary

The laws that govern natural systems can often be modelled mathematically using
partial differential equations (PDEs).
Usually the resultant PDEs are not amenable to analytical tools for
solving them, and numerical techniques remain the only resort in such cases.
As a result, efficiently solving PDEs numerically is a primary step towards
understanding the physics of various systems.
``SARAS`` is a general-purpose PDE solver written in object-oriented C++.
In ``SARAS``, the underlying mathematical constructs used to define a PDE, like
vector and scalar fields, are defined as classes.
Moreover, vector calculus operations associated with such fields, like gradient and
divergence, are defined as functions of the classes.

This design makes the code intuitive, allowing users to quickly cast PDEs into
readable code.
The initial conditions, boundary conditions, and source/forcing terms, which appear
commonly in many PDEs, are also defined as derived classes of base ``initial``,
``boundary`` and ``force`` classes.
These classes are written to be readily extensible, so that users can add custom
initial conditions, source terms, and so on.

``SARAS`` also includes solvers for hydrodynamic flow, namely the incompressible
Navier-Stokes initial value problem (IVP), as well as for scalar convection,
like Rayleigh Benard Convection.
Presently, we use semi-implicit Crank-Nicholson [@Crank:1947] method for time-advancing
the IVP.
The solver uses Marker and Cell method (MAC) [@Harlow:PF1965] for discretizing
the velocity and pressure fields.

# Mathematics

The Navier-Stokes equations, which govern the dynamics of fluid flow, can be written as
$$ \frac{\partial \mathbf{u}}{\partial t} + \mathbf{u}\cdot\nabla\mathbf{u} = -\nabla p + \mathbf{f} + \nu\nabla^2\mathbf{u}, $$
where $\mathbf{u}$ is the velocity field, $p$ is the pressure field, $\mathbf{f}$ is
the forcing term, and $\nu$ is the kinematic viscosity of the fluid.
The fluid is assumed incompressible. Hence $\nabla\cdot\mathbf{u} = 0$, and density is
assumed constant and equal to unity.

![For the flow simulation of decaying turbulence on a $256^3$ grid with ``TARANG`` (red lines) and ``SARAS`` (black-dashed lines):
  (a) plot of the total energy $E_u= \int d{\bf r} u^2/2$ vs $t$,
  (b) plot of $E_u(k)$ vs $k$ at $t =1$. \label{figure1}](figure1.png)

If the velocity and pressure field at time $t = n$ are denoted as $\mathbf{u}^n$ and $p^n$
respectively, then the corresponding fields at the next time-step, $t = n+1$, namely
$\mathbf{u}^{n+1}$ and $p^{n+1}$, can be calculated as described next [@Patankar:1972IJHMT].
Initially, an intermediate velocity field is calculated using the known values,
$\mathbf{u}^n$ and $p^n$, as
$$\mathbf{u}^* = \mathbf{u}_{n} + \Delta t\left[\nu\nabla^2 \left( \frac{\mathbf{u}_n + \mathbf{u}^*}{2}\right) - \mathbf{u}_n.\nabla\mathbf{u}_n - \nabla p_n\right]. $$
The forcing term has been neglected for simplicity.
Note that the diffusion term (also called the viscous term) has been split into two, with half
the contribution from the velocity field at $t = n$, and the other half attributed to the
guessed velocity field, $\mathbf{u}^*$.
This results in the following implicit equation,
$$\mathbf{u}^* - \Delta t\left[\frac{\nu\nabla^2\mathbf{u}^*}{2} \right ] = \mathbf{u}_{n} + \Delta t\left[\frac{\nu\nabla^2\mathbf{u}_n}{2} - \mathbf{u}_n.\nabla\mathbf{u}_n - \nabla p_n\right]. $$
The above equation has to be solved iteratively, and this is achieved through
OpenMP parallelized Jacobi iterations.
The intermediate velocity field, $\mathbf{u}^*$, will not satisfy the continuity equation,
and requires to be corrected appropriately.
This correction is obtained from the pressure correction term, which is in turn computed from
the pressure Poisson equation,
$$\nabla^2 p^* = \frac{\nabla.\mathbf{u}^*}{\Delta t}. $$
``SARAS`` uses a Geometric Multigrid library to solve the above equation [@Wesseling:MG2004; @Briggs:MG2000].
Presently the library offers the Full Multigrid (FMG) V-Cycle to solve the Poisson equation.
Other methods like F-Cycle and W-Cycle are planned updates to the library in future.

Finally, using this computed value of pressure correction, the velocity and pressure fields corresponding to the next time-step can now be obtained as
$$p^{n+1} = p^n + p^*,$$
$$\mathbf{u}^{n+1} = \mathbf{u}^* - \Delta t(\nabla p^*). $$
The solver also supports adaptive time-stepping, where the Courant-Friedrichs-Lewy (CFL) condition [@Courant:1928CFL] is used to dynamically compute the appropriate time-step
from the velocity field.

![For the flow simulation of decaying turbulence on a $256^3$ grid,
  vector plots of the velocity field and density plots of the vertical vorticity field ($\omega_z$) computed at the horizontal mid plane at $z=1/2$:
  for the data  from (a, c) ``TARANG``,  and (b, d) ``SARAS`` at $t=1$ and $t=3$. \label{figure2}](figure2.png)

# Results
We validate our code using two very well-known problems.
We simulate these problems using ``SARAS`` and compare the results with standard and validated solutions. 

## Problem 1
We simulate decaying turbulence using ``SARAS`` with Taylor-Green vortex [@Schranner:CP2013] as the initial condition.
The initial condition for the Taylor-Green velocity field is given by
$$
\mathbf{u}(x,y,z, t=0) = u_{0} \begin{bmatrix}
        \sin(2 \pi k_{0} x) \cos(2 \pi k_{0}y) \cos(2 \pi k_{0}z) \\
        -\cos(2 \pi k_{0} x) \sin(2 \pi k_{0}y) \cos(2 \pi k_{0}z) \\
        0
\end{bmatrix},
$$
where $u_0 = 1$ and  $k_0 = 1$.
We perform our simulation in a periodic box of size $1 \times 1 \times 1$ with a grid resolution of $256^{3}$, up to $t = 3.0$.
Here nondimensionalization of time is performed by $L/u_0$, where $L=1$.
The initial Reynolds number of the flow is $\mathrm{Re} = 1000$.
We choose a constant $dt = 0.001$ for time-integration.
Besides, we use a uniformly spaced mesh along all the three directions. 

Our results exhibit similarities with those obtained from the pseudo-spectral code, ``TARANG`` [@Chatterjee:JPDC2018].
The time-evolution of total kinetic energy ($\int d{\bf r} (u^2/2)$) is plotted in Figure \ref{figure1}(a).
The results from ``SARAS`` closely match the values from ``TARANG``,
with a maximum difference between the two energies of approximately $2.4\%$.
Besides, the flow profiles are also similar, as is evident from the density plots of Figure \ref{figure2}.
Here, the vertical component of vorticity, $\omega_z$, on the horizontal mid-plane is plotted at $t=1$ and $t=3$.
In Figure \ref{figure1}(b), we plot the energy spectrum at $t=1$.
The plot exhibits nearly similar multiscale evolution of the flow fields.
Interestingly, in both the results the energy spectrum in the inertial range is closer to the $k^{-5/3}$ prediction of Kolmogorov [@Kolmogorov:DANS1941Dissipation; @Kolmogorov:DANS1941Structure]. 

![Results from the simulation of lid-driven cavity on a $129^2$ grid with ``SARAS`` (orange lines), along with the data from [@Ghia:JCP1982] (blue stars):
  (a) The vertical profile of the x-component of velocity, $v_x$, along the line across the geometric center of the cavity
  (b) The horizontal profile of the y-compnent of velocity, $v_y$, along the line across the geometric center of the cavity. \label{figure3}](figure3.png)

## Problem 2
We solve the two-dimensional lid-driven cavity (LDC) problem using ``SARAS``,
and compare the results with those of [@Ghia:JCP1982].
LDC is an important fluid dynamic system serving as a benchmark for testing numerical methods.
The system comprises of a square cavity of dimension $1 \times 1$
with no-slip boundary conditions on all the four walls.
However, the top wall is moving laterally to the right with a constant velocity of $U = 1.0$,
and this serves as the reference velocity for non-dimensionalization of the problem.
The side of the cavity, $L = 1.0$, is the reference length.

At the start of the simulation, the velocity of the fluid is zero throughout the domain.
Thus, the fluid inside the cavity is driven impulsively by the top lid at the start of the simulation.
This results in the formation of a vortex at the upper-right corner of the cavity,
and this vortex rapidly grows in size and occupies the entire bulk of the cavity.
We simulate this on a $129 \times 129$ grid at $\mathrm{Re} = 1000$.
The result output by ``SARAS`` at $t = 30$ is used for comparison,
and the solver computes this solution with 4 MPI processes in under 12 minutes on a desktop workstation.
This case is also used as a test case for validation of the solver after its installation.

The bash and Python scripts used to run this test are supplied with the solver.
Upon running the test script, the solver is compiled and run using a
pre-defined set of parameters corresponding to the LDC problem.
After completion of the simulation, the shell script executes a Python script to read the output
from ``SARAS`` and compare the horizontal and vertical velocity profiles across
the geometric center of the square cavity.
The corresponding results from [@Ghia:JCP1982] are also available with the installation,
and this data is used by the Python script to show the match in profiles.

The result from this test is plotted in Figure \ref{figure3}.
We observe that the profile computed by ``SARAS`` matches the results from the literature very well.
This quick validation of the solver, which can be done as a part of its installation, is one of the strengths of this package.


# Acknowledgements

We acknowledge contributions from Gaurav Gautham, Saurav Bhattacharjee, and Rishabh Sahu,
and support from Prof Fahad Anwer during the genesis of this project.

---

# References

