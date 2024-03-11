
## Solving For Presure
### Solving Posssoins Equation
Solve Possoin's equation to find pressure, over all cells. Repeat until the solution converges 
(or until max interations have been reached).
$$
P_{i,j}^{k+1} = \frac{1}{4} \left( P_{i,j-1}^k + P_{i,j+1}^k + P_{i-1,j}^k + P_{i+1,j}^k - \frac{\rho \Delta x \Delta y}{\Delta t} \left( 2 \Delta x u_{i,j+1} - u_{i,j-1} + 2 \Delta y v_{i+1,j} - v_{i-1,j} \right) \right)
$$

By not mutating the current buffer whilst iterating over it (as the implicit SOR method does), 
concurrent solvers can be much easier and with requiring additional syncing.

### Boundary Conditions
Consider the **inlet** at $i=0$. Simple pressure gradient based in inlet velocity.
$$P_{0, j} = P_{1, j} - \frac{\rho V_{inlet} \Delta x} {2}$$

Consider the **Outlet** at $i=n$. Assumed boundary condition of 0 pressure.
$$P_{n, j} = 0$$

Consider a no-slip wall, copy adjacent cells pressure to force 0 pressure gradient. Velocity 
boundary conidtions are applied when evaulating the Navier Stokes equations to find the velocities.

## Finding Velocity
Once the pressure has beens solved at a given timestep, it can be applied to the Navier-Stokes equations to given the updated velocities for each cell.
$$
u_{i,j}^{k+1} = u_{i,j}^{k} - \frac{\Delta t}{\rho} \left( \frac{p_{i+1,j}^k - p_{i,j}^k}{\Delta x} \right) + \Delta t \nu \left( \frac{u_{i+1,j}^k - 2u_{i,j}^k + u_{i-1,j}^k}{\Delta x^2} + \frac{u_{i,j+1}^k - 2u_{i,j}^k + u_{i,j-1}^k}{\Delta y^2} \right) + \Delta t \cdot F_x
$$

$$
v_{i,j}^{k+1} = v_{i,j}^{k} - \frac{\Delta t}{\rho} \left( \frac{p_{i,j+1}^k - p_{i,j}^k}{\Delta y} \right) + \Delta t \nu \left( \frac{v_{i+1,j}^k - 2v_{i,j}^k + v_{i-1,j}^k}{\Delta x^2} + \frac{v_{i,j+1}^k - 2v_{i,j}^k + v_{i,j-1}^k}{\Delta y^2} \right) + \Delta t \cdot F_y
$$

Given that $\nu = \frac{\mu}{\rho}$ for incompressible flow.