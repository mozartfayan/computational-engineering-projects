# One-Dimensional Steady-State Convection-Diffusion Heat Transfer

## 1. Overview

This project presents a numerical solution of the one-dimensional steady-state convection-diffusion heat transfer problem using the Finite Volume Method (FVM).

Two spatial discretization schemes are implemented and compared:

- Central Difference Scheme (CDS)
- First-Order Upwind Scheme

The resulting linear systems are solved using the Tri-Diagonal Matrix Algorithm (TDMA). Numerical results are compared against the analytical solution.

## 2. Physical Problem

The problem considers steady-state heat transfer in a water flow with constant velocity through a one-dimensional domain.

**Physical properties and boundary conditions:**

| Parameter | Value |
|---|---:|
| Thermal conductivity (k) | 0.595 W/(m·K) |
| Density (ρ) | 997 kg/m³ |
| Domain length (L) | 0.1 m |
| Temperature at x = 0 | 150 °C |
| Temperature at x = L | 50 °C |

The governing equation is:

\[
k\frac{d^2T}{dx^2}-\rho u\frac{dT}{dx}=0
\]

## 3. Numerical Methodology

The computational domain is discretized into control volumes. The governing equation is integrated over each control volume, resulting in a system of algebraic equations.

Two convection discretization schemes are evaluated:

**Central Difference Scheme (CDS)**

The convective term is approximated using central differences. This scheme provides a more accurate representation under suitable mesh Peclet number conditions, but may produce non-physical oscillations when the mesh Peclet number is too high.

**First-Order Upwind Scheme**

The convective term is approximated using upstream values. This approach is generally more robust, although it introduces numerical diffusion.

**Linear System Solver**

Both discretization schemes produce tri-diagonal systems, which are solved using the Tri-Diagonal Matrix Algorithm (TDMA).

## 4. Analytical Solution

The numerical solutions are compared against the analytical temperature distribution:

\[
T(x)=T_0+(T_L-T_0)
\frac{e^{\rho ux/k}-1}{e^{\rho uL/k}-1}
\]

## 5. Computational Cases

Three simulation cases are considered:

| Case | Velocity (m/s) | Control Volumes |
|---|---:|---:|
| A | 0.01 | 5 |
| B | 0.10 | 5 |
| C | 0.10 | 20 |

These cases allow the investigation of the influence of flow velocity and mesh refinement on numerical accuracy and stability.

## 6. Results and Discussion

The numerical temperature profiles obtained using CDS and Upwind are compared with the analytical solution.

The mesh Peclet number is also evaluated to investigate the stability of the Central Difference Scheme and the influence of numerical diffusion introduced by the Upwind scheme.

### Temperature Distribution

<!-- Insert the temperature comparison plot here -->

## 7. Technologies

- Python
- NumPy
- Matplotlib
- Finite Volume Method
- TDMA

## 8. How to Run

Install the required Python libraries:

```bash
pip install numpy matplotlib
```

Run the simulation:

```bash
python main.py
```

## 9. Author

**Mozart Fayan Michel**

Environmental and Sanitary Engineer | PhD Researcher in Computational Fluid Dynamics (CFD)
