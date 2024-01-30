"""LU Decomposition

This module provides functionality to decompose given matrix in lower & upper triangular matrices via LU decomposition,

- lu_crout: LU decomposition via Crout's method.
- lu_dolittle LU decomposition via Dolittle's method.

Created on Jan 12 00:51:02 2024
"""

__all__ = ['lu_crout', 'lu_doolittle']

from custom_inherit import doc_inherit
from umatrix.matrix import Matrix, null_matrix

from .. import DOC_STYLE, IFloat, LMat, MatOrLList, N_DECIMAL
from ..__backend.matrix_ import round_matrix_


def lu_crout(matrix_a: MatOrLList) -> LMat:
    """Performs LU decomposition using Crout's method.

    Parameters
    ----------
    matrix_a:
        The matrix to decompose.


    Returns
    ------
        The LU decomposition following Crout's method.

    Raises
    ------
    NotASquareMatrix:
        If the matrix is not square.
    """

    matrix_a, lower_matrix, upper_matrix, n_rows = _lu_sanity_check(matrix_a)

    for i in range(n_rows):
        for j in range(i, n_rows):
            lower_sum = 0
            for k in range(i):
                lower_sum += lower_matrix[j][k] * upper_matrix[k][i]
            lower_matrix[j][i] = matrix_a[j][i] - lower_sum

        for j in range(i, n_rows):
            if i == j:
                upper_matrix[i][j] = 1
            else:
                upper_sum = 0
                for k in range(i):
                    upper_sum += lower_matrix[i][k] * upper_matrix[k][j]
                upper_matrix[i][j] = (matrix_a[i][j] - upper_sum) / lower_matrix[i][i]

    return [lower_matrix, upper_matrix]


@doc_inherit(lu_crout, style=DOC_STYLE)
def lu_doolittle(matrix_a: MatOrLList, n_decimal: IFloat = N_DECIMAL) -> LMat:
    """Performs LU decomposition using Doolittle's method.

    Returns
    ------
        The LU decomposition following Doolittle's method.
    """

    matrix_a, lower_matrix, upper_matrix, n_rows = _lu_sanity_check(matrix_a)

    for i in range(n_rows):
        upper_matrix[0][i] = matrix_a[0][i]

    for i in range(n_rows):
        lower_matrix[i][0] = matrix_a[i][0] / upper_matrix[0][0]

    for i in range(1, n_rows):
        for j in range(1, n_rows):
            if j == i:
                lower_matrix[i][j] = 1
            elif j < i:
                lower_sum = 0
                for k in range(j):
                    lower_sum += (lower_matrix[i][k] * upper_matrix[k][j])
                lower_matrix[i][j] = (matrix_a[i][j] - lower_sum) / upper_matrix[j][j]

        for j in range(1, n_rows):
            if j < i:
                upper_matrix[i][j] = 0
            else:
                upper_sum = 0
                for k in range(j):
                    upper_sum += (lower_matrix[i][k] * upper_matrix[k][j])
                upper_matrix[i][j] = matrix_a[i][j] - upper_sum

    lower_matrix = round_matrix_(lower_matrix, n_decimal)
    upper_matrix = round_matrix_(upper_matrix, n_decimal)

    return [lower_matrix, upper_matrix]


def _lu_sanity_check(matrix_a: Matrix):
    """
    Performs sanity check for LU decomposition.

    Parameters
    ----------
    matrix_a:
        The given matrix or list of list to be checked.

    Returns
    -------
        The given matrix, lower, and upper matrix, and number of rows in the matrix.
    """

    if not isinstance(matrix_a, Matrix):
        matrix_a = Matrix(matrix_a)

    n_rows, n_cols = matrix_a.n_rows, matrix_a.n_cols

    # if not matrix_a._is_symmetric():
    #     raise NotASquareMatrix('The given matrix is not square, LU decomposition cannot be performed.')

    lower_matrix = null_matrix(n_rows, n_cols)
    upper_matrix = null_matrix(n_rows, n_cols)

    return matrix_a, lower_matrix, upper_matrix, n_rows
