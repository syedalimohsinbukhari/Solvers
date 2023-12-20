"""Created on Dec 20 14:10:28 2023"""

# TODO: Add proper README.md


__all__ = ["LagrangeInterpolation",
           "newton_interpolation",
           "LinearSpline", "NaturalCubicSpline", "QuadraticSpline",
           "gauss_jacobi", "gauss_seidel",
           "newton_raphson_solver",
           "trapezoid_rule", 'simpson_13', 'simpson_38', 'boole', 'weddle_3h10', 'weddle_41h140',
           "euler_trapezoidal", "euler_method",
           "rk2_solver", 'rk3_solver', 'rk4_solver', 'rk4_multi_ode']

from .interpolators.LagrangeInterpolation import LagrangeInterpolation
from .interpolators.NewtonianInterpolators import newton_interpolation
from .interpolators.SplineInterpolation import LinearSpline, NaturalCubicSpline, QuadraticSpline
from .iterative_solvers.GaussJacobi import gauss_jacobi
from .iterative_solvers.GaussSeidel import gauss_seidel
from .iterative_solvers.NewtonCotesSolvers import (boole, newton_raphson_solver, simpson_13, simpson_38, trapezoid_rule,
                                                   weddle_3h10, weddle_41h140)
from .iterative_solvers.PredictorCorrectorMethods import euler_method, euler_trapezoidal
from .iterative_solvers.RungeKutta import rk2_solver, rk3_solver, rk4_multi_ode, rk4_solver
