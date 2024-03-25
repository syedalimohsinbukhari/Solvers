"""Cholesky Decomposition

This module provides the functionality to decompose a symmetric and positive definite matrix A using
Cholesky decomposition. The Cholesky decomposition expresses the matrix A as the product of a lower triangular matrix
L and its conjugate transpose L*.

Functions:
    - cholesky_decomposition: Performs Cholesky decomposition on the given matrix.
    - cholesky_ldl: Performs LDL Cholesky decomposition on the given matrix.

Created on Jan 09 02:03:16 2024
"""

__all__ = ['cholesky_decomposition', 'cholesky_ldl']

from math import sqrt

from custom_inherit import doc_inherit
from umatrix.matrix import Matrix, identity_matrix, null_matrix

from .. import DOC_STYLE, LList, LMat
from ..__backend.matrix_ import cholesky_sanity_check


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

    matrix_ = cholesky_sanity_check(matrix)

    n_rows, n_cols = matrix_.n_rows, matrix_.n_cols
    matrix_l = null_matrix(n_rows, n_cols)

    matrix_l[0][0] = sqrt(matrix_[0][0])

    for row_iter in range(1, n_rows):
        for col_iter in range(n_cols):
            if col_iter == 0:
                matrix_l[row_iter][col_iter] = matrix_[col_iter][row_iter] / matrix_l[0][0]
            elif col_iter < row_iter:
                temp_ = [matrix_l[row_iter][elem_] * matrix_l[col_iter][elem_] for elem_ in range(row_iter - 1)]
                matrix_l[row_iter][col_iter] = (matrix_[row_iter][col_iter] - sum(temp_)) / matrix_l[col_iter][col_iter]
            elif col_iter == row_iter:
                temp_ = sum((elem_**2 for elem_ in matrix_l[row_iter].elements))
                matrix_l[row_iter][col_iter] = sqrt(matrix_[row_iter][row_iter] - temp_)
            else:
                matrix_l[row_iter][col_iter] = 0

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

    def _ldl_reduction(working_index):
        # We don't bother about the row operations because the transformation is symmetric
        # This is done to reduce the number of operations to be performed on the matrix.
        for iter_ in range(working_index, n_rows):
            value_normalization = matrix_[iter_][working_index] / matrix_[working_index][working_index]

            matrix_l[working_index][iter_] = value_normalization if iter_ > 0 else matrix_[iter_][working_index]

            if iter_ > working_index:
                matrix_[iter_] -= value_normalization * matrix_[working_index]

        # extra step because we have to have 1 as the main diagonal element for matrix_l
        if working_index == 0:
            matrix_l[working_index][working_index] = 1

    matrix_ = cholesky_sanity_check(matrix)
    n_rows = matrix_.n_rows

    matrix_l = identity_matrix(n_rows)

    for row_number in range(n_rows):
        _ldl_reduction(row_number)

    return [matrix_l.t, matrix_.diagonal_of_matrix(), matrix_l]
