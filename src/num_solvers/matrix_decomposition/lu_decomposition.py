"""LU Decomposition

This module provides functionality to decompose given matrix in lower & upper triangular matrices via LU decomposition,

- lu_crout: LU decomposition via Crout's method.
- lu_dolittle LU decomposition via Dolittle's method.

Created on Jan 12 00:51:02 2024
"""

__all__ = ['lu_crout', 'lu_dolittle']

from custom_inherit import doc_inherit

from .cholesky_decomposition import cholesky_decomposition
from .. import DOC_STYLE, LLList, LList, N_DECIMAL
from ..__backend.errors_ import NotASquareMatrix
from ..__backend.matrix_decomposition_ import __lower_upper_matrices, __round_matrices


@doc_inherit(cholesky_decomposition, style=DOC_STYLE)
def lu_crout(matrix_a: LList, n_decimal: int = N_DECIMAL) -> LLList:
    """Performs LU decomposition using Crout's method.

    Returns
    ------
        The LU decomposition following Crout's method.

    Raises
    ------
    NotASquareMatrix:
        If the matrix is not square.
    """

    if len(matrix_a) != len(matrix_a[0]):
        raise NotASquareMatrix('The given matrix is not square, LU decomposition cannot be performed.')

    len_mat = len(matrix_a)

    lower_matrix, upper_matrix = __lower_upper_matrices(len_mat)

    for i in range(len_mat):

        for j in range(i, len_mat):
            lower_sum = 0
            for k in range(i):
                lower_sum += lower_matrix[j][k] * upper_matrix[k][i]
            lower_matrix[j][i] = matrix_a[j][i] - lower_sum

        for j in range(i, len_mat):
            if i == j:
                upper_matrix[i][j] = 1
            else:
                upper_sum = 0
                for k in range(i):
                    upper_sum += lower_matrix[i][k] * upper_matrix[k][j]
                upper_matrix[i][j] = (matrix_a[i][j] - upper_sum) / lower_matrix[i][i]

    return __round_matrices(lower_matrix, upper_matrix, n_decimal)


@doc_inherit(cholesky_decomposition, style=DOC_STYLE)
def lu_dolittle(matrix_a: LList, n_decimal: int = N_DECIMAL) -> LLList:
    """Performs LU decomposition using Dolittle's method.

    Returns
    ------
        The LU decomposition following Dolittle's method.

    Raises
    ------
    NotASquareMatrix:
        If the matrix is not square.
    """

    if len(matrix_a) != len(matrix_a[0]):
        raise NotASquareMatrix('The given matrix is not square, LU decomposition cannot be performed.')

    len_mat = len(matrix_a)

    lower_matrix, upper_matrix = __lower_upper_matrices(len_mat)

    for i in range(len_mat):
        upper_matrix[0][i] = matrix_a[0][i]

    for i in range(len_mat):
        lower_matrix[i][0] = matrix_a[i][0] / upper_matrix[0][0]

    for i in range(1, len_mat):
        for j in range(1, len_mat):
            if j == i:
                lower_matrix[i][j] = 1
            if j < i:
                lower_sum = 0
                for k in range(j):
                    lower_sum += (lower_matrix[i][k] * upper_matrix[k][j])
                lower_matrix[i][j] = (matrix_a[i][j] - lower_sum) / upper_matrix[j][j]

        for j in range(1, len_mat):
            if j < i:
                upper_matrix[i][j] = 0
            else:
                upper_sum = 0
                for k in range(j):
                    upper_sum += (lower_matrix[i][k] * upper_matrix[k][j])
                upper_matrix[i][j] = matrix_a[i][j] - upper_sum

    return __round_matrices(lower_matrix, upper_matrix, n_decimal)
