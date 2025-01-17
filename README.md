# 2D Cavity Flow with Explicit Projection Method: BiCGStab, FVM, Staggered Grid

This C++ code solves the incompressible Navier-Stokes equations for a two-dimensional cavity flow problem using the finite volume method.

## Overview

The code simulates the flow of a viscous, incompressible fluid within a closed cavity. It utilizes a time-stepping approach to solve the discretized Navier-Stokes equations. Key aspects of the implementation include:

* **Finite Volume Method:** The computational domain is divided into a grid of control volumes, and the governing equations are integrated over these volumes.

* **Staggered Grid:** Velocity components (u, v) and pressure (P) are stored at different locations within the control volume to avoid spurious oscillations.

* **Time Discretization:** An explicit time-stepping scheme is used for the momentum equations.

* **Explicit Projection Method:** The code employs an explicit projection method. First, a provisional velocity field is calculated, and then a pressure correction is applied to enforce incompressibility.

* **Boundary Conditions:** The code allows for the specification of various boundary conditions on the walls of the cavity, including velocity and temperature profiles.

* **Viscosity and Thermal Conductivity:** The code includes functions to calculate viscosity and thermal conductivity, although in the provided version, these are set to constant values.

* **Conjugate Gradient Solver:** The pressure Poisson equation is solved using the BiCGStab (Stabilized Bi-Conjugate Gradient) iterative method.

* **Temperature Equation:** The code includes the solution of the energy equation to simulate heat transfer within the cavity.

## Governing EquationsThe code solves the following equations:

* **Continuity Equation (Mass Conservation):** ∇ ⋅ **u** = 0

* **Momentum Equations (Newton's Second Law):** ∂**u**/∂t + (**u** ⋅ ∇) **u** = - (1/ρ) ∇P + ν ∇²**u**

* **Energy Equation (Heat Transfer):** ρ Cv (∂T/∂t + **u** ⋅ ∇T) = ∇ ⋅ (k ∇ T) + Φ where Φ represents the viscous dissipation term.

## Numerical Method

The code uses an explicit projection method to solve the Navier-Stokes equations:

1. **Provisional Velocity:** A provisional velocity field (**u***) is computed from the momentum equation without considering the pressure gradient:
    ∂**u**\*/∂t + (**u** ⋅ ∇) **u** = ν ∇²**u**
2. **Pressure Poisson Equation:** The divergence of the corrected velocity field is set to zero (incompressibility). This leads to a Poisson equation for the pressure:
    ∇²P = ρ/Δt  ∇ ⋅ **u***
3. **Pressure Boundary Condition:** The boundary condition for pressure is obtained by projecting the momentum equation onto the direction normal to the boundary and requiring that the corrected velocity equals the boundary velocity.
4. **Velocity Correction:** The provisional velocity is corrected using the calculated pressure gradient:
    **u**<sup>n+1</sup> = **u*** - (Δt/ρ) ∇P

## Code Structure
The code is organized into several functions, including:

* `main()`: The main function that controls the simulation flow.

* `initialization()`: Sets up the grid, allocates memory, and initializes variables.

* `calculate_boundary()`: Defines and applies boundary conditions.

* `provisional_velocity()`: Computes the intermediate velocity field (**u***).

* `pressure_iteration()`: Solves the pressure Poisson equation.

* `new_velocity()`: Updates the velocity field based on the pressure correction.

* `temperature_update()`: Solves the energy equation.

* `write_results()`: Outputs the simulation results to files.

* `calculate_shear()`: Calculate shear stress and vorticity.

* `interpolate_boundary()`: Interpolates boundary conditions for a given side and location.

## Compilation

To compile the code, you will need a C++ compiler (like g++). Navigate to the directory containing HS_cavityVar.cpp in your terminal and use the following command:
```bash
g++ HS_cavityVar.cpp -o HS_cavityVar -lm
```

**Note:**

* `mem_alloc.cpp` is assumed to be a separate file containing memory allocation/deallocation functions used by the code. Ensure this file is in the same directory or adjust the compilation command accordingly.
* `-lm` links the math library, which is necessary for functions like `sqrt`, `sin`, `pow`, etc.

## Running the Simulation

After successful compilation, you can run the executable using the following command:
```bash
./HS_cavityVar
```

The simulation will proceed, and output files containing the results will be generated in the same directory.

## Output Files

The code generates several output files:

* `Value_UV.dat`: Contains the x and y coordinates, u and v velocities, stress, vorticity, pressure, temperature, viscosity, and thermal conductivity at each cell center.
* `Left_wall_stress.dat`: Stores data related to the left wall, including shear stress components, temperature, viscosity, thermal conductivity, and the top plate velocity profile.
* `Right_wall_stress.dat`: Stores data related to the right wall.
* `Wall_uw.dat`: Contains wall velocity gradients and related quantities.
* `Gen_info.txt`: Provides general simulation parameters and calculated values like shear force.
* `recordLeft.dat`, `recordBottom.dat`, `recordRight.dat`, `recordTop.dat`: These files likely contain the applied boundary conditions.

## Parameters

The code defines several parameters at the beginning of the `HS_cavityVar.cpp` file that control the simulation:

* `N`: Number of grid cells in each direction.
* `MAX_ERR`: Maximum error tolerance for the pressure solver.
* `EPS`: A small value used for initialization.
* `ERR_UV`: Error tolerance for velocity convergence.
* `PR_ITE`: Maximum number of iterations for the pressure solver.
* `L`: Physical length of the cavity.
* `RHO`: Density of the fluid.
* `DELT`: Time step size.
* `ITE`: Total number of time iterations.
* `temperature`: Initial temperature.
* `Sigma`, `AtomMass`, `CV`: Parameters related to fluid properties.

You can modify these parameters to explore different flow regimes and simulation setups.

## Memory Allocation

The code uses custom memory allocation functions (`allocate` and `deallocate`) defined in `mem_alloc.h` and implemented in `mem_alloc.cpp`. These functions are used for dynamic allocation of 2D and 3D arrays.

## Further Information

For a deeper understanding of the numerical methods employed, you may want to research:

* Finite Volume Method for Fluid Dynamics
* BiCGStab Method
* Explicit Projection Method
* Discretization schemes for convection and diffusion terms

This README provides a basic overview of the `HS_cavityVar.cpp` code. Understanding the underlying principles of computational fluid dynamics will be beneficial for interpreting the results and modifying the code.