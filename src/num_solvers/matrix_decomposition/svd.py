"""Created on Jan 18 00:48:38 2024"""

from operator import mul

from . import LMat
from .matrix import Matrix, identity_matrix, null_matrix, vector_mag
from .qr_decomposition import qr_decomposition
from .. import FList, Func, IFloat, N_DECIMAL, OptIFloat, TOLERANCE
from ..__backend.extra_ import filter_similar_values, linear_list


# TODO: See if numpy can be removed

# taken from https://www-users.cse.umn.edu/~olver/aims_/qr.pdf
def eigen_values_qr(matrix: Matrix, n_iter: int = 1000) -> LMat:
    # if matrix.is_singular:
    #     raise ValueError('Cannot use QR decomposition method for eigen values on singular matrices.')

    id_ = identity_matrix(matrix.n_rows, matrix.n_cols)

    for _ in range(n_iter):
        temp_ = qr_decomposition(matrix)
        id_ *= temp_[0]
        matrix = mul(*temp_[::-1])

    return [id_, matrix.diagonal()]


def rayleigh_quotient_method(matrix: Matrix, n_iter: int = 50, initial_guess: int = 200) -> list:

    identity_ = identity_matrix(matrix.n_rows, matrix.n_cols)

    b_matrix = null_matrix(n_iter, 1).elements
    mu_matrix = null_matrix(n_iter, 1).elements

    b_matrix[0], mu_matrix[0] = Matrix([[1]] * matrix.n_rows), initial_guess

    for _ in range(n_iter - 1):
        f1 = (matrix - mu_matrix[_] * identity_).inverse() * b_matrix[_]
        b_matrix[_ + 1] = f1 / vector_mag(f1)

        mu_matrix[_ + 1] = b_matrix[_ + 1].t * matrix * b_matrix[_ + 1]
        mu_matrix[_ + 1] /= b_matrix[_ + 1].t * b_matrix[_ + 1]

    return mu_matrix[-1]


def eigen_value_multi(iter_func: Func, matrix: Matrix, range_: tuple = (-10, 10),
                      n_iter: int = 50) -> list:
    l_list = linear_list(range_[0], range_[1], 2 * sum(map(abs, range_)))

    eigen_value = [iter_func(matrix, n_iter, mu_value) for mu_value in l_list]

    return list(sorted(filter_similar_values(eigen_value, TOLERANCE), reverse=True))


def eigen_identity(n_rows: int, diag_value, n_cols: OptIFloat = None) -> Matrix:
    null_ = null_matrix(n_rows, n_cols if n_cols else n_rows)

    for i in range(null_.n_rows):
        null_[i][i] = diag_value

    return null_


def characteristic_polynomial(matrix: Matrix) -> FList:
    return faddeev_le_verrier(matrix)


def faddeev_le_verrier(matrix: Matrix) -> FList:
    new_matrix_ = [0] * (matrix.n_cols + 1)
    coefficient = [0] * (matrix.n_cols + 1)

    new_matrix_[0], coefficient[0] = identity_matrix(matrix.n_rows), 1

    for i in range(1, matrix.n_rows + 1):
        temp_ = matrix * new_matrix_[i - 1]
        coefficient[i] = (-1 / i) * sum(temp_.diagonal().elements)
        new_matrix_[i] = temp_ + identity_matrix(matrix.n_rows, value=coefficient[i])

    return coefficient


def svd(matrix: Matrix, n_decimal: IFloat = N_DECIMAL):
    raise NotImplementedError()
