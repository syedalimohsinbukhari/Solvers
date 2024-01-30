"""QR decomposition, Gram-Schmidt

This module provides functionality to perform QR decomposition of the given matrix following Gram-Schmidt method via,

- qr_decomposition

Created on Jan 12 12:42:02 2024
"""

__all__ = ['qr_decomposition']

from umatrix.matrix import Matrix, null_matrix, vector_mag

from .. import LMat, MatOrLList
from ..__backend.matrix_ import remove_zeroed_columns


def qr_decomposition(matrix: MatOrLList) -> LMat:
    """
    Performs the QR decomposition for the given matrix.

    Parameters
    ----------
    matrix:
        The given matrix to perform QR decomposition on.

    Returns
    -------
        Q and R decomposed matrices according to Gram-Schmidt method.
    """

    if matrix.n_cols == 1:
        raise ValueError("Can not perform QR decomposition for mx1 matrix.")

    u_matrix = Matrix([[i[j] for i in matrix.elements] for j in range(matrix.n_cols)])

    v_matrix = null_matrix(1, matrix.n_cols).elements
    q_matrix = null_matrix(1, matrix.n_cols).elements
    r_matrix = null_matrix(matrix.n_rows, matrix.n_cols)

    v_matrix[0] = u_matrix[0]

    temp_ = []
    for elements in range(1, matrix.n_cols):
        for k in range(elements):
            dot_ = u_matrix[elements].dot(v_matrix[k]) / vector_mag(v_matrix[k], True)
            temp_.append(dot_ * v_matrix[k])
        v_matrix[elements] = u_matrix[elements] - sum(temp_)
        temp_ = []

    # v_matrix = [round_matrix_(i, n_decimal) for i in v_matrix]

    for index, value in enumerate(v_matrix):
        if all(x == 0 for x in value):
            q_matrix[index] = Matrix([0] * len(value))
        else:
            q_matrix[index] = value / vector_mag(value)

    q_matrix = Matrix([i.elements for i in q_matrix]).t

    for row in range(matrix.n_rows):
        for col in range(matrix.n_cols):
            if row > col:
                r_matrix[row][col] = 0
            else:
                r_matrix[row][col] = u_matrix[col].dot(q_matrix.t[row])

    q_matrix = remove_zeroed_columns(q_matrix.t).t
    r_matrix = remove_zeroed_columns(r_matrix)

    return [q_matrix, r_matrix]
