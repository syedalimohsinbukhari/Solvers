"""Custom errors module

This module contains custom error classes for interpolation methods. The error classes include,
- InterpolationException: The base exception class for interpolation methods.
- AtLeastOneParameterRequired: This error is raised when neither ``function`` not ``function_values`` parameter is
  supplied in an interpolation class.
- WrongBoundaryEquation: This error is raised when the user provides wrong string for boundary equation.

Created on Jan 07 16:05:26 2024
"""

__all__ = ['AtLeastOneParameterRequired', 'WrongBoundaryEquation']


class InterpolationException(BaseException):
    pass


class AtLeastOneParameterRequired(InterpolationException):
    pass


class WrongBoundaryEquation(InterpolationException):
    pass
