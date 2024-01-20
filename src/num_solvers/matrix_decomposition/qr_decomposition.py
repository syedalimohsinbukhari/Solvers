"""QR decomposition, Gram-Schmidt

This module provides functionality to perform QR decomposition of the given matrix following Gram-Schmidt method via,

- qr_decomposition

Created on Jan 12 12:42:02 2024
"""

__all__ = ['qr_decomposition']

from . import LMat, N_DECIMAL
from .matrix import Matrix, null_matrix, vector_mag
from ..__backend.extra_ import round_matrix_


def schur_decomposition(matrix: Matrix, n_decimal: int = N_DECIMAL) -> LMat:
    raise NotImplementedError()


def qr_decomposition(matrix: Matrix, n_decimal: int = N_DECIMAL) -> LMat:
    """
    Performs the QR decomposition for the given matrix.

    Parameters
    ----------
    matrix:
        The given matrix to perform QR decomposition on.
    n_decimal:
        Number of digits to round off to. Default is 8.

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

    v_matrix = [round_matrix_(i, n_decimal) for i in v_matrix]

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

    q_matrix = _remove_zero_columns(q_matrix.t).t
    r_matrix = _remove_zero_columns(r_matrix)

    q_matrix = round_matrix_(q_matrix, n_decimal)
    r_matrix = round_matrix_(r_matrix, n_decimal)

    return [q_matrix, r_matrix]


def _remove_zero_columns(matrix_to_modify: Matrix) -> Matrix:
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
            f1 = new_[:-cond]
            return Matrix(f1)
    else:
        return Matrix(new_)
