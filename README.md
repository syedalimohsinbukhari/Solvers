# num_solvers

---
A library for some basic numerical solvers.

## Why another numerical solver library?

---
The library is the product of my PhD coursework for **Advance Numerical Techniques** by _**Dr. Umair**_ at _Institute of Space Technology, Pakistan, Islamabad_ and is not intended to replace any existing numerical solving packages. This is for educational purpose only.

Most of the modules included in this library have been implemented efficiently in the mainstream packages like **NumPy**, **SciPy**, and more. The users are encouraged to use those packages for their daily tasks.

## What's in this library?

---
This library includes implementations of,
* [Newton-Cotes integrators](src/num_solvers/integrators/newton_cotes_integrators.py),
  * Trapezoid Rule
  * Simpson's 1/3 Rule
  * Simpson's 3/8 Rule
  * Boole's Method
  * Weddle 3/10 Rule
  * Weddle 4/140 Rule