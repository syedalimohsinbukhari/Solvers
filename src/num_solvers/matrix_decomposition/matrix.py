"""Created on Oct 04 21:36:44 2023"""

from umatrix.matrix import Matrix, identity_matrix

from .. import FList


# TODO: Check the setting of values inside column matrices, they're acting up


def characteristic_polynomial(matrix: Matrix) -> FList:
    return faddeev_le_verrier(matrix)


def faddeev_le_verrier(matrix: Matrix) -> FList:
    new_matrix_: list = [0] * (matrix.n_cols + 1)
    coefficient: list = [0] * (matrix.n_cols + 1)

    new_matrix_[0], coefficient[0] = identity_matrix(matrix.n_rows), 1

    for i in range(1, matrix.n_rows + 1):
        temp_ = matrix * new_matrix_[i - 1]
        coefficient[i] = (-1 / i) * sum(temp_.diagonal().elements)
        new_matrix_[i] = temp_ + identity_matrix(matrix.n_rows, value=coefficient[i])

    return coefficient
