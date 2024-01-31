"""
Cholesky Decomposition

This module provides the functionality to decompose a symmetric and positive definite matrix A using
Cholesky decomposition. The Cholesky decomposition expresses the matrix A as the product of a lower triangular matrix
L and its conjugate transpose L*.

Functions:
    - cholesky_decomposition: Performs Cholesky decomposition on the given matrix.

Created on Jan 09 02:03:16 2024
"""

__all__ = ['cholesky_decomposition']

from math import sqrt

from umatrix.matrix import Matrix, null_matrix

from .. import LList, LMat
from ..__backend.errors_ import NonSymmetricMatrix, NotPositiveDefinite
from ..__backend.matrix_ import reduce_to_zeros


def cholesky_decomposition(matrix: Matrix or LList) -> LMat:
    """
    Performs the Cholesky decomposition on the given matrix.

    Parameters
    ----------
    matrix:
        The matrix to decompose.

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

    if not isinstance(matrix, Matrix):
        matrix = Matrix(matrix)

    if not matrix.is_symmetric():
        raise NonSymmetricMatrix('The matrix is not symmetric. Can not perform Cholesky decomposition.')

    if not matrix.is_positive_definite():
        raise NotPositiveDefinite('The matrix is not positive definite. Can not perform Cholesky decomposition.')

    n_rows, n_cols = matrix.n_rows, matrix.n_cols
    matrix_l = null_matrix(n_rows, n_cols)

    matrix_l[0][0] = sqrt(matrix[0][0])

    for i in range(1, n_rows):
        for k in range(n_rows):
            if k == 0:
                matrix_l[i][k] = matrix[k][i] / matrix_l[0][0]
            elif k < i:
                f1_ = [matrix_l[i][j] * matrix_l[k][j] for j in range(i - 1)]
                matrix_l[i][k] = (matrix[i][k] - sum(f1_)) / matrix_l[k][k]
            elif k == i:
                f2_ = sum((elem**2 for elem in matrix_l[i].elements))
                matrix_l[i][k] = sqrt(matrix[i][i] - f2_)
            else:
                matrix_l[i][k] = 0

    matrix_l = reduce_to_zeros(matrix_l)

    return [matrix_l, matrix_l.t]
