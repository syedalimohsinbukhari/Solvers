"""Steepest descent

This module provides functionality to solve given linear system of equations using steepest descent method via,

- steepest_descent :math:`^{[1]}`
- modified_steepest_descent :math:`^{[2]}`

Notes
-----
Prefer **modified_steepest_descent** over **steepest_descent** for quick convergence.

References
----------
:math:`^{[1]}` https://www.phys.uconn.edu/~rozman/Courses/m3511_18s/downloads/steepest-descent.pdf

:math:`^{[2]}` https://doi.org/10.1186/s13662-020-02715-9

Created on Jan 22 23:00:22 2024
"""

from custom_inherit import doc_inherit

from .matrix import Matrix, vector_mag
from .. import DOC_STYLE, IFloat, N_DECIMAL, TOLERANCE
from ..__backend.extra_ import round_matrix_


# taken from https://www.phys.uconn.edu/~rozman/Courses/m3511_18s/downloads/steepest-descent.pdf
def steepest_descent(matrix: Matrix, solution: Matrix, initial_guess: Matrix, n_iterations: int = 250,
                     tolerance: IFloat = TOLERANCE, n_decimal: int = N_DECIMAL) -> Matrix:
    """
    Solve a linear system using the steepest descent method.

    Parameters
    ----------
    matrix : Matrix
        Coefficient matrix of the linear system.
    solution : Matrix
        Right-hand side vector of the linear system.
    initial_guess : Matrix
        Initial guess for the solution vector.
    n_iterations : int, optional
        Maximum number of iterations. Default is 50.
    tolerance : IFloat, optional
        Tolerance for convergence. Default is TOLERANCE.
    n_decimal : int, optional
        Number of decimal places to round the result. Default is N_DECIMAL.

    Returns
    -------
    Matrix:
        The solution vector.

    Notes
    -----
        * This function uses the steepest descent method to iteratively solve the linear system :math:`Ax = b`.
        * For larger matrices, it is advised to increase the number of iterations to get proper results.

    References
    ----------
    https://www.phys.uconn.edu/~rozman/Courses/m3511_18s/downloads/steepest-descent.pdf
    """

    r: list = [0] * n_iterations
    new_guess: list = [0] * n_iterations
    new_guess[0] = initial_guess

    processed_iters = 0

    for iter_ in range(1, 25):
        r[iter_ - 1] = solution - (matrix * new_guess[iter_ - 1])
        alpha = (r[iter_ - 1].t * r[iter_ - 1]) / (r[iter_ - 1].t * matrix * r[iter_ - 1])
        new_guess[iter_] = new_guess[iter_ - 1] + alpha * r[iter_ - 1]

        temp_ = vector_mag(new_guess[iter_] - new_guess[iter_ - 1])

        if temp_ < tolerance:
            processed_iters = iter_
            break

    return round_matrix_(new_guess[processed_iters], n_decimal)


# direct link: https://d-nb.info/1215094116/34
@doc_inherit(steepest_descent, style=DOC_STYLE)
def modified_steepest_descent(matrix: Matrix, solution: Matrix, initial_guess: Matrix, n_iterations: int = 250,
                              tolerance: IFloat = TOLERANCE, n_decimal: int = N_DECIMAL):
    """Solve a linear system using the modified steepest descent method.

    References
    ----------
    "The steepest descent of gradient-based iterative method for solving rectangular linear systems with an \
application to Poisson’s equation." https://doi.org/10.1186/s13662-020-02715-9
    """

    def iter_summation(calculation_type: str = '1'):
        if calculation_type == '1':
            condition = matrix.n_cols
            mat1, mat2 = c_, c
        else:
            condition = matrix.n_rows
            mat1, mat2 = d_, d

        sum1, sum2, cond2 = 0, 0, matrix.n_cols

        for i in range(condition):
            sum2 += (mat2[i] - sum(mat1[i][j] * new_guess[iter_ - 1][j] for j in range(matrix.n_cols))).elements[0]**2

        return sum2

    tau: list = [0] * n_iterations

    new_guess: list = [0] * n_iterations
    new_guess[0] = initial_guess

    c, c_ = matrix.t * solution, matrix.t * matrix
    d, d_ = matrix * c, matrix * c_

    processed_iters = 0

    for iter_ in range(1, n_iterations):
        tau[iter_ - 1] = iter_summation() / iter_summation('2')
        new_guess[iter_] = new_guess[iter_ - 1] + tau[iter_ - 1] * (c - c_ * new_guess[iter_ - 1])

        temp_ = vector_mag(new_guess[iter_] - new_guess[iter_ - 1])

        if temp_ < tolerance:
            processed_iters = iter_
            break

    return round_matrix_(new_guess[processed_iters], n_decimal)
