"""Custom errors module

This module contains custom error classes for interpolation methods. The error classes include,
- InterpolationException: The base exception class for interpolation methods.
- AtLeastOneParameterRequired: This error is raised when neither ``function`` not ``function_values`` parameter is
  supplied in an interpolation class.
- WrongBoundaryEquation: This error is raised when the user provides wrong string for boundary equation.

Created on Jan 07 16:05:26 2024
"""

__all__ = ['AtLeastOneParameterRequired', 'WrongBoundaryEquation', 'NonSymmetricMatrix', 'NotPositiveDefinite']


class SolversErrors(BaseException):
    pass


class InterpolationException(SolversErrors):
    pass


class AtLeastOneParameterRequired(SolversErrors):
    pass


class WrongBoundaryEquation(SolversErrors):
    pass


class NonSymmetricMatrix(SolversErrors):
    pass


class NotPositiveDefinite(SolversErrors):
    pass
