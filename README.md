# Project for Computing at Scale

--- 

## Visualizing an Ordinary Differential Equation with perturbed term

---

*This software aims to simplify the task of solving ODE's with the help of Perturbation Theory. It also helps visualizing the effect of the perturbed term in an equation.*

### The need for this Software

This software helps Non-mathematical background people to apply perturbation theory on Linear and Non linear ODEs to understand the effects of different terms in a system. 

### Inputs and Outputs of the software

- **Inputs**: The input would be a non-dimensional equation with an option whether the equation is linear or non linear. Later they also have to input the value of the perturbed parameter e and the time range or number of steps. The equation would be in the form of ax''+ bx'+ ec = 0. The user has to put in values for a,b and c initially followed by e and number of time steps, tfinal. 

- **Outputs**: The output of the system would be a graph for typical parameter and saves a file with numerical solution of the equation for the number of steps given.

---

### Useful Libraries

This software will use libraries like odeint and rk4_solver for solving equations. ODEint is a package that helps solve differential equations numerically. It used RK4 method to solve the equations. Documentation for the same can be found here: https://headmyshoulder.github.io/odeint-v2/examples.html.  Another way of solving this is using rk4() library. Documentation for it can be found here: https://people.math.sc.edu/Burkardt/cpp_src/rk4/rk4.html. 
