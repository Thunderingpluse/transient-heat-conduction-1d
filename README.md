# 1D Transient Heat Conduction Solver

## Aim
To solve the 1D transient heat conduction equation in a solid slab using the explicit Finite Difference Method (FTCS scheme) with an insulated boundary and a constant temperature boundary.

## Theory
The governing equation for 1D transient heat conduction without heat generation is:

$$\frac{\partial T}{\partial t} = \alpha \frac{\partial^2 T}{\partial x^2}$$

Where $\alpha$ is the thermal diffusivity ($m^2/s$).

### Discretization (FTCS Explicit Scheme)
Using a forward difference in time and central difference in space:

$$\frac{T_i^{n+1} - T_i^n}{\Delta t} = \alpha \frac{T_{i+1}^n - 2T_i^n + T_{i-1}^n}{\Delta x^2} \implies T_i^{n+1} = \beta T_{i+1}^n + (1 - 2\beta) T_i^n + \beta T_{i-1}^n$$

Where $\beta = \frac{\alpha \Delta t}{\Delta x^2}$ is the grid Fourier number. 
For **numerical stability**, the Fourier number must satisfy:
$$\beta \le 0.5$$

### Boundary Conditions
- **Left Wall ($x=0$)**: Insulated boundary (adiabatic, $\frac{\partial T}{\partial x} = 0$). Using a second-order ghost node substitution gives:
  $$T_0^{n+1} = 2\beta T_1^n + (1 - 2\beta) T_0^n$$
- **Right Wall ($x=L$)**: Prescribed constant temperature (Dirichlet, $T_{N-1}^{n+1} = T_{boundary}$).

## File Structure
- `transient.py` - The explicit iterative solver simulating heat diffusion over time steps and plotting individual node temperatures.
- `output.txt` - Step-by-step logs showing temperatures across the nodes for successive time-steps until convergence.
- `graph1.png` - Visual layout containing grid plots showing temperature variation versus time-step for each node.

## How to Run
Ensure you have the required dependencies:
```bash
pip install numpy matplotlib
python transient.py
```
