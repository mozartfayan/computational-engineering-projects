# Transient 1D Heat Conduction: Finite Volume Method and TDMA

## Overview

This repository contains a Python implementation for solving a one-dimensional transient heat conduction problem in a solid bar using the **Finite Volume Method (FVM)** and the **Tri-Diagonal Matrix Algorithm (TDMA)**.

The main objective is to numerically solve the transient heat diffusion equation and evaluate the accuracy of the numerical solution by comparing it with its corresponding analytical solution, obtained through a Fourier series expansion.

This project was developed as part of a Computational Fluid Dynamics (CFD) numerical methods exercise, with an emphasis on numerical discretization, iterative time integration, linear system solving, and error analysis.

## Physical Problem

The model considers transient heat conduction along a one-dimensional solid bar with the following boundary conditions:

- **Left boundary (x = 0):** Adiabatic boundary condition (zero heat flux).
- **Right boundary (x = L):** Prescribed temperature (Dirichlet boundary condition).
- **Initial condition:** Uniform temperature distribution throughout the bar.

The governing equation is the one-dimensional transient heat conduction equation:

\[
\rho c_p \frac{\partial T}{\partial t}
=
k \frac{\partial^2 T}{\partial x^2}
\]

where:

- \(T\) is temperature;
- \(\rho\) is density;
- \(c_p\) is specific heat capacity;
- \(k\) is thermal conductivity;
- \(t\) is time;
- \(x\) is the spatial coordinate.

## Numerical Methodology

The numerical implementation consists of three main components:

1. **Finite Volume Method (FVM):** Spatial discretization of the governing heat conduction equation using control volumes.
2. **Fully Implicit Time Integration:** Temporal discretization using a fully implicit formulation.
3. **TDMA Solver:** Solution of the resulting tridiagonal linear system at each time step.

The numerical solution is computed at every time step, while the analytical solution is evaluated at the same spatial locations and time instants.

## Analytical Solution

The analytical solution is obtained using a Fourier series expansion, considering the same initial and boundary conditions as the numerical model.

This solution serves as a reference for evaluating the numerical method's accuracy.

## Error Analysis

The code calculates the absolute error between the numerical and analytical temperature distributions at each time step.

Two error metrics are evaluated:

- **Maximum Absolute Error:** Maximum absolute temperature difference across all control volumes.
- **Mean Absolute Error:** Average absolute temperature difference across all control volumes.

The results are stored throughout the simulation, allowing the temporal evolution of numerical accuracy to be investigated.

## Results and Visualization

The script automatically generates two figures:

1. **Temperature Distribution:** Comparison between numerical and analytical temperature profiles at selected time instants.
2. **Absolute Error Evolution:** Maximum and mean absolute errors as functions of time.

All generated figures are saved in the `out/` directory.

## Technologies and Libraries

- Python
- NumPy
- Matplotlib

## How to Run

### 1. Clone the repository

```bash
git clone https://github.com/YOUR_USERNAME/YOUR_REPOSITORY.git
cd YOUR_REPOSITORY
```

### 2. Install dependencies

```bash
python -m pip install numpy matplotlib
```

### 3. Run the simulation

```bash
python solucao_tdma_lista6_versao3.py
```

The numerical results will be printed in the terminal, and the generated figures will be saved in the `out/` directory.

## Project Objectives

This project aims to:

- Implement the Finite Volume Method for a transient heat conduction problem.
- Develop and apply the TDMA algorithm to solve tridiagonal systems.
- Compare numerical and analytical solutions.
- Investigate numerical accuracy through absolute error metrics.
- Develop practical skills in scientific computing and numerical methods for engineering applications.

## Author

**Mozart Fayan Michel**

Environmental and Sanitary Engineer | PhD Researcher in Computational Fluid Dynamics (CFD)

**Research interests:** Computational Fluid Dynamics, Numerical Methods, Heat Transfer, Non-Newtonian Fluid Mechanics, and Scientific Computing.
