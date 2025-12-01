# num_solvers

[![Python](https://img.shields.io/badge/python-3.9%20%7C%203.10%20%7C%203.11-blue)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A Python library for basic numerical solvers, providing implementations of various numerical methods for integration, interpolation, root-finding, and matrix decomposition.

## Table of Contents

- [About the Project](#about-the-project)
- [Features](#features)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Modules Overview](#modules-overview)
- [Use Cases](#use-cases)
- [Documentation](#documentation)
- [Requirements](#requirements)
- [License](#license)

## About the Project

**num_solvers** is an educational numerical methods library developed as part of PhD coursework for **Advanced Numerical Techniques** by _**Dr. Umair**_ at the _Institute of Space Technology, Pakistan, Islamabad_.

### Why another numerical solver library?

This library is **not intended to replace** existing numerical computing packages. It is designed for:

- **Educational purposes**: Learn and understand the implementation details of classical numerical methods
- **Reference implementations**: Clear, readable code that demonstrates how these algorithms work
- **Academic use**: Suitable for coursework, research, and teaching numerical methods

> **Note**: For production use, you are encouraged to use well-established packages like **NumPy**, **SciPy**, and similar libraries that offer optimized and thoroughly tested implementations.

## Features

- **Numerical Integration**: Newton-Cotes integrators including Trapezoidal, Simpson's, Boole's, and Weddle's rules
- **Interpolation Methods**: Newton, Lagrange, and Spline interpolation techniques
- **Root Finding**: Multiple algorithms for finding roots of functions and polynomials
- **ODE Solvers**: Runge-Kutta methods (RK2, RK3, RK4) for ordinary differential equations
- **Matrix Decomposition**: LU, QR, Cholesky, SVD, and other decomposition methods
- **System Solvers**: Gauss-Jacobi, Gauss-Seidel, and predictor-corrector methods

## Installation

### Using pip

```bash
pip install num_solvers
```

### From source

```bash
# Clone the repository
git clone https://github.com/syedalimohsinbukhari/Solvers.git
cd Solvers

# Install in development mode
pip install -e .
```

## Quick Start

Here are some examples to get you started:

### Numerical Integration

```python
from num_solvers.integrators import trapezoidal_rule, simpson_one_third

# Define a function to integrate
def f(x):
    return x**2

# Integrate using the trapezoidal rule
result = trapezoidal_rule(f, a=0, b=1, n=100)
print(f"Trapezoidal rule result: {result}")

# Integrate using Simpson's 1/3 rule
result = simpson_one_third(f, a=0, b=1, n=100)
print(f"Simpson's 1/3 rule result: {result}")
```

### Interpolation

```python
from num_solvers.interpolators import lagrange_interpolation

# Define data points
x_data = [0, 1, 2, 3]
y_data = [1, 2, 0, 5]

# Interpolate at a specific point
x_interp = 1.5
y_interp = lagrange_interpolation(x_data, y_data, x_interp)
print(f"Interpolated value at x={x_interp}: {y_interp}")
```

### Root Finding

```python
from num_solvers.iterative_solvers import bisection_method, newton_raphson

# Define a function and its derivative
def f(x):
    return x**3 - 2*x - 5

def f_prime(x):
    return 3*x**2 - 2

# Find root using bisection method
root = bisection_method(f, a=2, b=3, tol=1e-6)
print(f"Root (bisection): {root}")

# Find root using Newton-Raphson method
root = newton_raphson(f, f_prime, x0=2.5, tol=1e-6)
print(f"Root (Newton-Raphson): {root}")
```

## Modules Overview

### 1. Integrators

Numerical integration methods for computing definite integrals.

| Method | Description |
|--------|-------------|
| Trapezoidal Rule | Basic numerical integration using trapezoidal approximation |
| Simpson's 1/3 Rule | Integration using quadratic polynomial approximation |
| Simpson's 3/8 Rule | Integration using cubic polynomial approximation |
| Boole's Rule | Higher-order Newton-Cotes formula |
| Weddle's Rule | 6-point Newton-Cotes formula (3/10 and 41/140 variants) |

📖 [Detailed documentation](docs/integrators/newton-cotes-integrators.md)

### 2. Interpolation Methods

Methods for estimating values between known data points.

| Method | Description |
|--------|-------------|
| Newton Forward Interpolation | For equally spaced data near the beginning |
| Newton Backward Interpolation | For equally spaced data near the end |
| Newton Divided Difference | For unequally spaced data |
| Lagrange Interpolation | Polynomial interpolation through all points |
| Linear Spline | Piecewise linear interpolation |
| Quadratic Spline | Piecewise quadratic interpolation |
| Natural Cubic Spline | Smooth cubic interpolation |

📖 Detailed documentation: [Newton](docs/interpolators/newton-interpolators.md) | [Lagrange](docs/interpolators/lagrange-interpolators.md) | [Spline](docs/interpolators/spline-interpolators.md)

### 3. Iterative Solvers

Methods for solving equations and systems of equations.

#### ODE Solvers (Runge-Kutta Methods)

| Method | Description |
|--------|-------------|
| RK2 (Heun's Method) | Second-order Runge-Kutta method |
| RK3 | Third-order Runge-Kutta method |
| RK4 | Classical fourth-order Runge-Kutta method |

📖 [Detailed documentation](docs/iterative-solvers/rk-solvers.md)

#### Root Finding Methods

| Method | Description |
|--------|-------------|
| Bisection Method | Bracketing method for continuous functions |
| False Position (Regula Falsi) | Improved bracketing method |
| Newton-Raphson | Fast convergence using derivatives |
| Secant Method | Approximates Newton's method without derivatives |
| Muller's Method | Uses quadratic interpolation |
| Ridder's Method | Exponential interpolation-based bracketing |
| Sidi's Method | Generalized secant method |
| Steffensen's Method | Self-accelerating method |
| Laguerre's Method | Specialized for polynomial roots |

📖 [Detailed documentation](docs/iterative-solvers/root-finding--functions.md)

#### System of Equations Solvers

| Method | Description |
|--------|-------------|
| Gauss-Jacobi | Iterative method for diagonally dominant systems |
| Gauss-Seidel | Improved iterative method with faster convergence |
| Predictor-Corrector | Multi-step method for ODEs |

### 4. Matrix Decomposition Methods

Methods for decomposing matrices into simpler forms.

| Method | Description |
|--------|-------------|
| Cholesky Decomposition | For symmetric positive-definite matrices |
| LU Decomposition (Crout) | Lower-Upper decomposition using Crout's method |
| LU Decomposition (Doolittle) | Lower-Upper decomposition using Doolittle's method |
| QR Decomposition | Orthogonal-triangular decomposition |
| Householder Transformation | For QR decomposition and reduction |
| Givens Rotation | For QR decomposition via plane rotations |
| Singular Value Decomposition (SVD) | General matrix factorization |
| Gauss Elimination | Direct method for solving linear systems |
| Steepest Descent | Iterative optimization method |
| Rayleigh Quotient | For eigenvalue problems |

📖 Detailed documentation available in the [docs/matrix-decomposition/](docs/matrix-decomposition/) directory.

## Use Cases

### 1. Educational Learning
- Understand how numerical algorithms work under the hood
- Compare different methods for the same problem
- Study convergence properties and error analysis

### 2. Scientific Computing
- Numerical integration for physics simulations
- Solving ODEs for dynamical systems
- Interpolation for data analysis and curve fitting

### 3. Engineering Applications
- Root-finding for nonlinear equations in circuit analysis
- Matrix decomposition for structural analysis
- Numerical integration for signal processing

### 4. Academic Research
- Reference implementations for algorithm comparison
- Baseline methods for developing new algorithms
- Teaching numerical methods courses

## Documentation

Detailed documentation for each module is available in the [docs/](docs/) directory:

- **Integrators**: [Newton-Cotes Integrators](docs/integrators/newton-cotes-integrators.md)
- **Interpolators**: 
  - [Newton Interpolators](docs/interpolators/newton-interpolators.md)
  - [Lagrange Interpolators](docs/interpolators/lagrange-interpolators.md)
  - [Spline Interpolators](docs/interpolators/spline-interpolators.md)
- **Iterative Solvers**:
  - [Runge-Kutta Solvers](docs/iterative-solvers/rk-solvers.md)
  - [Root Finding (Functions)](docs/iterative-solvers/root-finding--functions.md)
  - [Root Finding (Polynomials)](docs/iterative-solvers/root-finding--polynomials.md)
- **Matrix Decomposition**: See [docs/matrix-decomposition/](docs/matrix-decomposition/)

## Requirements

- Python >= 3.9
- NumPy < 2.0.1
- setuptools
- custom-inherit
- umatrix

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Author

**Syed Ali Mohsin Bukhari**  
📧 syedali.b@outlook.com

---

⭐ If you find this library useful for educational purposes, consider giving it a star!
