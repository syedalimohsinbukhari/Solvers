"""Created on Jan 28 15:43:08 2024"""

from umatrix.matrix import Matrix, identity_matrix, vector_mag

from .. import N_DECIMAL
from ..__backend.matrix_ import reduce_to_zeros, round_matrix_


def givens_rotation(matrix: Matrix, n_decimal: int = N_DECIMAL, override_original: bool = False):
    matrix_ = matrix if override_original else Matrix(matrix.elements[:])

    givens_ = identity_matrix(matrix_.n_rows)

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
        givens_ *= identity_.t

    givens_ = round_matrix_(reduce_to_zeros(givens_), n_decimal)
    matrix_ = round_matrix_(reduce_to_zeros(matrix_), n_decimal)

    return givens_, matrix_
