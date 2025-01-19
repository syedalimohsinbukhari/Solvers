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

# __all__ = ['AtLeastOneParameterRequired', 'WrongBoundaryEquation', 'NonSymmetricMatrix', 'NotASquareMatrix',
#            'NumStepCanNotBeFloat', 'NotPositiveDefinite', 'DegreeOfPolynomialNotCorrect', 'IndexCanNotBeSlice',
#            'InconsistentDimensions', 'XFindNotDefined']


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


#######################################################################################################################

class DivisionByMatrix(NumSolversErrors):
    def __init__(self, message="Can't divide a matrix by matrix. Invalid Operation."):
        super().__init__(message)


class MatrixDimensionsMismatch(NumSolversErrors):
    def __init__(self, message=''):
        m = "Matrix dimensions do not match for multiplication.\n"
        super().__init__(m + message)


class NotASquareMatrixError(NumSolversErrors):
    def __init__(self, message="The provided matrix is not a square matrix.\nCan't perform LU Decomposition."):
        super().__init__(message)


class SlicingNotAllowed(NumSolversErrors):
    def __init__(self, message="Slicing a matrix is not possible, please use integer indices."):
        super().__init__(message)


class IndexOutOfBounds(NumSolversErrors):
    def __init__(self, message="The given index doesn't exist for the matrix."):
        super().__init__(message)


class DeterminantIsZero(NumSolversErrors):
    def __init__(self, message="The given matrix is singular and its inverse can't be calculated."):
        super().__init__(message)
