"""Matrix_

This module provides a helpful functionality to display the matrices in fractions via,

- :class:`InFractions`

Created on Jan 10 00:01:13 2024
"""

__all__ = ['map_to_matrix', 'remove_zeroed_columns', 'round_matrix_', 'reduce_to_zeros']

from umatrix.matrix import Matrix

from .core_helpers_ import round_list_
from .. import IFloat, N_DECIMAL, TOLERANCE


def reduce_to_zeros(matrix: Matrix, tolerance: IFloat = TOLERANCE) -> Matrix:
    """
    Replaces the values in the matrices below tolerance to 0s.

    Parameters
    ----------
    matrix:
        The matrix in which the values are to be replaced with zeros.
    tolerance:
        The tolerance level below which the values are reduced to zero. Default is 1e-8.

    Returns
    -------
        The matrix with 0s instead of very small values.
    """

    for i, row in enumerate(matrix.elements):
        for j, element in enumerate(row):
            if abs(element) < tolerance:
                matrix[i][j] = 0

    return matrix


def remove_zeroed_columns(matrix_to_modify: Matrix) -> Matrix:
    """
    Removes the extra columns from the final Q and R decomposed matrices.

    Parameters
    ----------
    matrix_to_modify:
        The decomposed Q or R matrix.

    Returns
    -------
        Matrix with removed columns after all 0s.
    """

    new_, all_zeros = [0] * matrix_to_modify.n_rows, []

    for i, elem in enumerate(matrix_to_modify.elements):
        if not all(j == 0 for j in elem):
            new_[i] = elem
        else:
            new_[i] = elem
            all_zeros.append(i)

    if bool(all_zeros):
        if len(new_) == all_zeros[0]:
            return Matrix(new_)
        else:
            cond = len(matrix_to_modify) - all_zeros[0]
            return Matrix(new_[:-cond])
    else:
        return Matrix(new_)


def round_matrix_(matrix: Matrix, n_decimal: int = N_DECIMAL) -> Matrix:
    """
    Maps the round function to a matrix.

    Parameters
    ----------
    matrix:
        Given matrix to round off.
    n_decimal:
        Number of decimal places to round off to.

    Returns
    -------
        Rounded off matrix.
    """

    return Matrix(round_list_(matrix.elements, n_decimal))


def map_to_matrix(function, matrix: Matrix):
    for i in range(matrix.n_rows):
        for j in range(matrix.n_cols):
            matrix[j][j] = function(matrix[i][j])

    return matrix
