"""Cholesky decomposition

This module provides the functionality to decompose a symmetric and positive definite matrix A, using
Cholesky decomposition.

Created on Jan 09 02:03:16 2024
"""

from math import sqrt

from .. import LLList, LList
from ..__backend.errors_ import NonSymmetricMatrix, NotPositiveDefinite

__all__ = ['cholesky_decomposition']


def cholesky_decomposition(matrix_a: LList, n_decimal: int = 4) -> LLList:
    """
    Performs the Cholesky decomposition on the given matrix.

    Parameters
    ----------
    matrix_a:
        The matrix to decompose.
    n_decimal:
        The number of digits to round off for the result.

    Returns
    -------
        The Cholesky decomposition of the matrix. Matrix L, and L_star.

    Raises
    ------
    NonSymmetricMatrix:
        If the matrix is not symmetric.
    NotPositiveDefinite:
        If the matrix is not positive definite.
    """

    if not __is_symmetric(matrix_a):
        raise NonSymmetricMatrix('The matrix is not symmetric. Can not perform Cholesky decomposition.')

    if not __is_positive_definite(matrix_a):
        raise NotPositiveDefinite('The matrix is not positive definite. Can not perform Cholesky decomposition.')

    matrix_l = [[0] * len(matrix_a) for _ in range(len(matrix_a))]

    n_dimensions = len(matrix_a)

    matrix_l[0][0] = sqrt(matrix_a[0][0])
    for i in range(1, n_dimensions):
        for k in range(n_dimensions):
            if k == 0:
                matrix_l[i][k] = matrix_a[k][i] / matrix_l[0][0]
            elif k < i:
                f1 = matrix_a[i][k]
                f1_ = [matrix_l[i][j] * matrix_l[k][j] for j in range(i - 1)]
                matrix_l[i][k] = (f1 - sum(f1_)) / matrix_l[k][k]
            elif k == i:
                f4 = sum((i**2 for i in matrix_l[i]))
                matrix_l[i][k] = sqrt(matrix_a[i][i] - f4)
            else:
                matrix_l[i][k] = 0

    matrix_l = [[round(element, n_decimal) for element in row] for row in matrix_l]
    matrix_l_star = [[matrix_l[row][col] for row in range(n_dimensions)] for col in range(n_dimensions)]

    return [matrix_l, matrix_l_star]


def __is_symmetric(matrix: LList) -> bool:
    """Checks if the given matrix is symmetric or not.

    Parameters
    ----------
    matrix:
        The matrix to check.

    Returns
    -------
        True if the matrix is symmetric, else False.
    """

    n_dim = len(matrix)
    matrix_transpose = [[matrix[row][col] for row in range(n_dim)] for col in range(n_dim)]

    return matrix == matrix_transpose


def __is_positive_definite(matrix: LList) -> bool:
    """Checks if the matrix is positive definite or not.

    Parameters
    ----------
    matrix:
        The matrix to check.

    Returns
    -------
        True if the matrix is positive definite, else False.

    """

    n_dimensions = len(matrix)

    # Check if the matrix is square
    if len(matrix) != len(matrix[0]):
        return False

    # Check if the matrix is symmetric
    if not __is_symmetric(matrix):
        return False

    # Check if all leading principal minors have positive determinants
    for i in range(1, n_dimensions + 1):
        minor = [row[:i] for row in matrix[:i]]
        determinant = sum((minor[j][j] for j in range(i)))
        if determinant <= 0:
            return False

    return True
