# Project for Computing at Scale

--- 

## Solving a system of stiff, coupled, non-linear equations using IVP solver for coupled oscillators.  

---

*This code aims to simplify the task of solving ODE's numerically and understand the system of equations better. It also helps visualizing the effects of each term in a system of equations.*

### The need for this code

This code is used in solving system of equation for a non-linear electro-mechanical coupled oscillator. This oscillator mimics a part of an low powered sensor. This sensor is understood to improve the energy transfer from electrical part to the mechanical part of the system. The system of equations are as follows:
---
                                                    my" + ky - 1/2 q^2/eA = kg
                                                    Lq" + q/C0 + qy/eA = 0
---
The Non-dimensional version of this is:
---
                                                    y" + y + 1/4 C1q^2 = 1/2
                                                    q" + C2q + 2C3qy = 0
---                                                    
Here C1, C2 and C3 are constant parameters which define the system's internal relations. The goal is the toggle these parameters and the boundary conditions [y0 ,q0 ,y0',q0'] to see the nature of the system. These parameters have internal relation - C2+C3 = 1 to equate natural frequencies of both the oscillators. 

---

### Inputs and Outputs of the software

- **Inputs**: The input would be two non-dimensional parameters [C1, C2] and boundary conditions. It would also ask for number of time steps that this system needs to be running.
  The format of input should be
  < C1 > < C2 > <time_steps> <[y0 ,q0 ,y0',q0']>

- **Outputs**: The output of the system would be a graphs for x(t), q(t), x'(t) and q'(t) and saves a file with numerical solution of the equation for the number of steps given.

---

### Useful Libraries

This software will use libraries like odeint and rk4_solver for solving equations. ODEint is a package that helps solve differential equations numerically. It used RK4 method to solve the equations. Documentation for the same can be found here: https://headmyshoulder.github.io/odeint-v2/examples.html.  Another way of solving this is using rk4() library. Documentation for it can be found here: https://people.math.sc.edu/Burkardt/cpp_src/rk4/rk4.html. Sundials is also a good option to solve the ODE's numerically. 
