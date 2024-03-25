"""Givens rotation

Provides functionality to apply Givens rotation to the matrix, via,

- given_rotation

Created on Jan 28 15:43:08 2024
"""

from umatrix.matrix import Matrix, identity_matrix, matrix_copy, vector_mag

from .. import LMat


def givens_rotation(matrix: Matrix, overwrite: bool = False) -> LMat:
    """
    Apply Givens rotations to triangularize a matrix.

    Parameters
    ----------
    matrix:
        The input matrix to be triangularized.
    overwrite:
        If True, the original matrix is modified in place. If False (default),
        a new matrix is created, and the original matrix remains unchanged.

    Returns
    -------
    LMat
        A list containing two matrices:

            1. The Givens rotation matrix (product of all Givens rotations).
            2. The triangularized matrix.
    """

    matrix_ = matrix_copy(matrix, overwrite)
    givens_matrix = identity_matrix(matrix_.n_rows)

    for i in range(matrix_.n_rows - 1):
        identity_ = identity_matrix(matrix_.n_rows)

        r_magnitude = vector_mag(matrix_.t[i].elements[i:i + 2])
        cosine = matrix_[i][i] / r_magnitude
        sine = -matrix[i + 1][i] / r_magnitude

        identity_[i][i] = cosine
        identity_[i + 1][i + 1] = cosine
        identity_[i][i + 1] = -sine
        identity_[i + 1][i] = sine

        matrix_ = identity_ * matrix_
        givens_matrix *= identity_.t

    return [givens_matrix, matrix_]
