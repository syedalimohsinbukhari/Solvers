"""Custom errors module

This module contains custom error classes for interpolation methods. The error classes include,

- NumSolversErrors:
    The base exception class for this module.

# Interpolation Errors
- AtLeastOneParameterRequired:
    When neither ``function`` not ``function_values`` parameter is supplied for interpolation.
- WrongBoundaryEquation:
    When the user provides wrong string for boundary equation.

# Matrix Errors
- NonSymmetricMatrix:
    When the matrix provided isn't symmetric.
- NotPositiveDefinite:
    When the matrix provided isn't positive definite.

# Polynomial Errors
- DegreeOfPolynomialNotCorrect:
    When length of the polynomial and the degree doesn't match.

# NumSolversMatrix Errors
- IndexCanNotBeSlice:
    When the index to be set is slice instead of integer.
- InconsistentDimensions:
    When the members of the array have inconsistent length.

Created on Jan 07 16:05:26 2024
"""

__all__ = ['AtLeastOneParameterRequired', 'WrongBoundaryEquation', 'NonSymmetricMatrix', 'NotASquareMatrix',
           'NumStepCanNotBeFloat', 'NotPositiveDefinite', 'DegreeOfPolynomialNotCorrect', 'IndexCanNotBeSlice',
           'InconsistentDimensions', 'XFindNotDefined']


class NumSolversErrors(BaseException):
    pass


class AtLeastOneParameterRequired(NumSolversErrors):
    pass


class WrongBoundaryEquation(NumSolversErrors):
    pass


class NumStepCanNotBeFloat(NumSolversErrors):
    pass


class XFindNotDefined(NumSolversErrors):
    pass


class NonSymmetricMatrix(NumSolversErrors):
    pass


class NotASquareMatrix(NumSolversErrors):
    pass


class NotPositiveDefinite(NumSolversErrors):
    pass


class DegreeOfPolynomialNotCorrect(NumSolversErrors):
    pass


class IndexCanNotBeSlice(NumSolversErrors):
    pass


class InconsistentDimensions(NumSolversErrors):
    pass
