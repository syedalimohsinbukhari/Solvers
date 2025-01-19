"""Created on Feb 13 00:31:34 2024"""

__all__ = ['gauss_jacobi', 'gauss_seidel']

from custom_inherit import doc_inherit
from umatrix.matrix import Matrix, identity_matrix, vector_mag

from .. import DOC_STYLE, FList, IFloat, MAX_ITER, TOLERANCE
from ..__backend.matrix_ import gauss_methods_sanity_check, upper_diagonal_lower_matrices, upper_lower_matrices


def gauss_seidel(matrix: Matrix, solution: FList, initial_guess: FList, max_iterations: int = MAX_ITER,
                 tolerance: int = TOLERANCE) -> Matrix:
    """
    Solves a system of linear equations using the Gauss-Seidel iterative method.

    Parameters
    ----------
    matrix:
        Coefficient matrix of the system.
    solution:
        Vector representing the constant terms in the equations.
    initial_guess:
        Initial guess for the solution vector.
    max_iterations:
        Maximum number of iterations. Default is 1000.
    tolerance:
        Tolerance for convergence. Default is 1e-8.

    Returns
    -------
    Matrix
        Solution vector.

    Notes
    -----
    1. The method works for diagonally dominant matrices properly.
       For a 3x3 matrix, the diagonally dominant matrix is defined as

       :math:`|a_{11}| \\geq |a_{12}| + |a_{13}|`

       :math:`|a_{22}| \\geq |a_{21}| + |a_{23}|`

       :math:`|a_{33}| \\geq |a_{31}| + |a_{32}|`

    2. This method also works for symmetric and positive definite matrices.
    """

    matrix_, solution, initial_guess = gauss_methods_sanity_check(matrix,
                                                                  solution,
                                                                  initial_guess)

    matrix_l, matrix_u = upper_lower_matrices(matrix_)

    initial_guess: list = [initial_guess]

    for i in range(max_iterations):
        matrix_t = -matrix_l.inverse() * matrix_u
        matrix_c = matrix_l.inverse() * solution

        initial_guess.append(matrix_t * initial_guess[-1] + matrix_c)

        temp_ = vector_mag(matrix_ * initial_guess[-1] - solution)
        if temp_ < tolerance:
            break

    return initial_guess[-1]


@doc_inherit(gauss_seidel, style=DOC_STYLE)
def gauss_jacobi(matrix: Matrix, solution: FList, initial_guess: FList, max_iterations: int = MAX_ITER,
                 tolerance: int = TOLERANCE) -> Matrix:
    """Solves a system of linear equations using the Gauss-Jacobi iterative method."""

    matrix_, solution, initial_guess = gauss_methods_sanity_check(matrix,
                                                                  solution,
                                                                  initial_guess)

    matrix_l, matrix_d, matrix_u = upper_diagonal_lower_matrices(matrix_)

    initial_guess: list = [initial_guess]

    for i in range(max_iterations):
        matrix_t = -matrix_d.inverse() * (matrix_l + matrix_u)
        matrix_c = matrix_d.inverse() * solution

        initial_guess.append(matrix_t * initial_guess[-1] + matrix_c)

        temp_ = vector_mag(matrix_ * initial_guess[-1] - solution)
        if temp_ < tolerance:
            break

    return initial_guess[-1]


@doc_inherit(gauss_seidel, style=DOC_STYLE)
def weighted_gauss_jacobi(matrix: Matrix, solution: FList, initial_guess: FList, max_iterations: int = MAX_ITER,
                          tolerance: int = TOLERANCE, weight: IFloat = 2 / 3) -> Matrix:
    """Solves a system of linear equations using the weighted Gauss-Jacobi iterative method."""

    matrix_, solution, initial_guess = gauss_methods_sanity_check(matrix,
                                                                  solution,
                                                                  initial_guess)

    matrix_l, matrix_d, matrix_u = upper_diagonal_lower_matrices(matrix_)

    weight, identity_ = weight, identity_matrix(matrix_.n_rows)

    initial_guess: list = [initial_guess]

    for i in range(max_iterations):
        matrix_t = weight * matrix_d.inverse() * solution

        f1_ = weight * matrix_d.inverse() * matrix_
        matrix_c = (identity_ - f1_) * initial_guess[-1]

        initial_guess.append(matrix_t + matrix_c)

        temp_ = vector_mag(matrix_ * initial_guess[-1] - solution)
        if temp_ < tolerance:
            break

    return initial_guess[-1]
