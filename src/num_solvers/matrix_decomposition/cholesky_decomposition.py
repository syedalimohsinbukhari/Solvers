"""Cholesky decomposition

This module provides the functionality to decompose a symmetric and positive definite matrix A, using
Cholesky decomposition.

- cholesky_decomposition: Performs matrix decomposition using Cholesky method.

Created on Jan 09 02:03:16 2024
"""

__all__ = ['cholesky_decomposition']

from math import sqrt

from .. import LLList, LList, N_DECIMAL
from ..__backend.errors_ import NonSymmetricMatrix, NotPositiveDefinite
from ..__backend.extra_ import round_list_
from ..__backend.matrix_decomposition_ import __is_positive_definite, __is_symmetric, __transpose


def cholesky_decomposition(matrix_a: LList, n_decimal: int = N_DECIMAL) -> LLList:
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

    n_dimensions = len(matrix_a)
    matrix_l = [[0] * n_dimensions for _ in range(n_dimensions)]

    matrix_l[0][0] = sqrt(matrix_a[0][0])

    for i in range(1, n_dimensions):
        for k in range(n_dimensions):
            if k == 0:
                matrix_l[i][k] = matrix_a[k][i] / matrix_l[0][0]
            elif k < i:
                f1_ = [matrix_l[i][j] * matrix_l[k][j] for j in range(i - 1)]
                matrix_l[i][k] = (matrix_a[i][k] - sum(f1_)) / matrix_l[k][k]
            elif k == i:
                f2_ = sum((i**2 for i in matrix_l[i]))
                matrix_l[i][k] = sqrt(matrix_a[i][i] - f2_)
            else:
                matrix_l[i][k] = 0

    matrix_l = [round_list_(row, n_decimal) for row in matrix_l]
    matrix_l_star = __transpose(matrix_l)

    return [matrix_l, matrix_l_star]
