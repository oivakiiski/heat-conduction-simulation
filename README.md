# Heat Conduction Simulation in MATLAB

This project was completed as part of my Applied Mathematics 1 course at LUT University.

The objective was to numerically model transient heat conduction through an
axisymmetric cylindrical aluminium domain and determine the resulting
temperature distribution and heat flux.

## Method

The model discretises the cylindrical domain in the radial and axial
directions and solves the heat equation numerically using the
Crank–Nicolson method.

The implementation includes:

- Temperature-dependent physical problem represented using aluminium
  material properties
- 2D axisymmetric spatial discretisation
- Sparse Laplacian operator
- Crank–Nicolson time integration
- Fixed-temperature boundary conditions
- Iterative convergence to steady state
- Heat-flux calculation using Fourier's law
- Visualisation of temperature distribution and heat-flux vectors

## Results

The simulation produces a two-dimensional temperature distribution through
the cylindrical domain and calculates the corresponding radial and axial
heat-flux components.

## What I learned

This project developed my understanding of translating a physical engineering
problem into a numerical model. It required combining heat-transfer theory,
numerical methods, matrix operations and MATLAB programming to create and
evaluate a working simulation.

## Why it matters

Although the physical system is different from vehicle performance, the
project demonstrates my ability to model a physical system, implement a
numerical simulation in MATLAB, analyse the resulting data and evaluate
system behaviour.
