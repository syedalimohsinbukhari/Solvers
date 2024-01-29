"""Created on Jan 25 06:38:05 2024"""

from umatrix.matrix import Matrix, identity_matrix, null_matrix, vector_mag


# taken from https://www-users.cse.umn.edu/~olver/aims_/qr.pdf


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
