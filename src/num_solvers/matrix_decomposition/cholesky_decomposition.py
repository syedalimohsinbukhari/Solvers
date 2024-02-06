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

import copy
from math import sqrt

from custom_inherit import doc_inherit
from umatrix.matrix import Matrix, identity_matrix, null_matrix

from .. import DOC_STYLE, LList, LMat
from ..__backend.matrix_ import cholesky_sanity_check, reduce_to_zeros


def cholesky_decomposition(matrix: Matrix or LList) -> LMat:
    """
    Performs the Cholesky decomposition on the given matrix.

    Parameters
    ----------
    matrix:
        The matrix to decompose.

    Returns
    -------
    LMat:
        The Cholesky decomposition of the matrix. Matrix L, and L_star.
    """

    matrix = cholesky_sanity_check(matrix)

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


# taken from https://services.math.duke.edu/~jdr/2021f-218/materials/week11.pdf
@doc_inherit(cholesky_decomposition, style=DOC_STYLE)
def cholesky_ldl(matrix: Matrix or LList) -> LMat:
    """

    Returns
    -------
    LMat:
        The LDL Cholesky decomposition of the matrix. Matrix L, diagonal matrix D, and L_star.
    """

    def _ldl_reduction(r_index):
        for iter_ in range(r_index, matrix_.n_rows):
            value_normalization = matrix_[iter_][r_index] / matrix_[r_index][r_index]

            matrix_l[r_index][iter_] = value_normalization if iter_ > 0 else matrix_[iter_][r_index]

            if iter_ > r_index:
                matrix_[iter_] -= value_normalization * matrix_[r_index]

        # extra step
        if r_index == 0:
            matrix_l[0][0] = 1

    matrix = cholesky_sanity_check(matrix)

    matrix_ = Matrix(copy.deepcopy(matrix.elements[:]))
    matrix_l = identity_matrix(matrix_.n_rows)

    for row_number in range(matrix_.n_rows):
        _ldl_reduction(row_number)

    return [matrix_l.t, matrix_.diagonal_of_matrix(), matrix_l]
