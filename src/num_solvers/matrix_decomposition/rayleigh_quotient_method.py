"""Rayleigh Quotient

This module provides an implementation of the Rayleigh Quotient method for finding the eigenvalue
of a matrix. The Rayleigh Quotient is a scalar value that characterizes the eigenvalues of a matrix.

References: https://www-users.cse.umn.edu/~olver/aims_/qr.pdf

Created on Jan 25 06:38:05 2024
"""

from umatrix.matrix import Matrix, identity_matrix, null_matrix, vector_mag

from .. import IFloat, MatOrLList


def rayleigh_quotient_method(matrix: MatOrLList, n_iter: int = 50, initial_guess: int = 200) -> IFloat:
    """
    Apply the Rayleigh Quotient method to estimate the dominant eigenvalue of a matrix.

    Parameters
    ----------
    matrix : MatOrLList
        The input matrix for which the eigenvalue is to be estimated.
    n_iter : int, optional
        The number of iterations for the Rayleigh Quotient method. Default is 50.
    initial_guess : int, optional
        The initial guess for the dominant eigenvalue. Default is 200.

    Returns
    -------
    IFloat
        The estimated dominant eigenvalue.

    Notes
    -----
    The Rayleigh Quotient method iteratively refines the estimate of the dominant eigenvalue
    by updating the eigenvector and eigenvalue.

    References
    ----------
    https://www-users.cse.umn.edu/~olver/aims_/qr.pdf
    """

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
