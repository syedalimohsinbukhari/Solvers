"""Created on Dec 20 14:10:28 2023"""

# TODO: Add proper README.md


__all__ = ['NewtonCotesIntegrators',
           'LagrangeInterpolation', 'NewtonianInterpolationMethods', 'SplineInterpolation',
           'GaussSeidel', 'GaussJacobi', 'PredictorCorrectorMethods', 'RungeKutta',
           'polynomial', 'non_polynomial']

from .integrators import NewtonCotesIntegrators
from .interpolators import LagrangeInterpolation, NewtonianInterpolationMethods, SplineInterpolation
from .iterative_solvers import non_polynomial, polynomial
from .iterative_solvers.odes import RungeKutta
from .iterative_solvers.system_of_equations import GaussJacobi, GaussSeidel, PredictorCorrectorMethods
