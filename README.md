# num_solvers

A library for some basic numerical solvers.

## Why another numerical solver library?

The library is the product of my PhD coursework for **Advance Numerical Techniques** by _**Dr. Umair**_ at _Institute of
Space Technology, Pakistan, Islamabad_ and is not intended to replace any existing numerical solving packages. This is
for educational purpose only.

Most of the modules included in this library have been implemented efficiently in the mainstream packages like **NumPy
**, **SciPy**, and more. The users are encouraged to use those packages for their daily tasks.

## What's in this library?

This library includes implementations of,

1. Integrators
    * [Newton-Cotes integrators](https://github.com/syedalimohsinbukhari/Solvers/blob/md-docs/docs/newton-cotes-integrators.md)
        * Trapezoidal rule
        * Simpson 1/3 rule
        * Simpson 3/8 rule
        * Boole's rule
        * Weddle 3/10 rule
        * Weddle 41/140 rule
2. Interpolation methods
    * [Newton interpolation method](https://github.com/syedalimohsinbukhari/Solvers/blob/md-docs/docs/newton-interpolators.md)
        * Forward interpolation
        * Backward interpolation
        * Divided difference interpolation
    * [Lagrange's interpolation method](https://github.com/syedalimohsinbukhari/Solvers/blob/md-docs/docs/lagrange-interpolators.md)
    * [Spline interpolation](https://github.com/syedalimohsinbukhari/Solvers/blob/md-docs/docs/spline-interpolators.md)
        * Linear spline
        * Quadratic spline
        * Natural cubic spline
3. Iterative solvers
    * ODE solvers
        * [Runge-Kutta solvers](https://github.com/syedalimohsinbukhari/Solvers/blob/md-docs/docs/rk-solvers.md)
    * [System of equation solvers]()
        * Gauss-Jacobi solver
        * Gauss-Seidel solver
        * Predictor-corrector solver
    * [Function root]()
        * Bisection method
        * False position method
        * Generalized secant method
        * Muller method
        * Newton-Raphson method
        * Ridder method
        * Sidi method
        * Steffensen method
    * [Polynomial root]()
        * Laguerre method
4. Matrix decomposition methods
    * [Cholesky decomposition]()
    * [Extra]()
    * [Gauss elimination]()
    * [Given's rotation]()
    * [Householder method]()
        * Householder reduction
        * Householder QR decomposition
    * [LU decomposition]()
        * LU Crout
        * LU Doolittle
    * [QR decomposition]()
    * [Rayleigh quotient method]()
    * [Steepest descent method]()
        * Steepest descent
        * Modified steepest descent
    * [Singular value decomposition]()