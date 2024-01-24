"""Created on Jan 21 15:05:20 2024"""

from functools import reduce
from math import sqrt
from operator import mul

from . import LMat, Matrix
from .matrix import identity_matrix, vector_mag
from .. import IFloat, TOLERANCE


def householder_reduction(matrix: Matrix, overwrite_original: bool = False) -> LMat:
    """
    Reduces the given matrix using Householder's method.

    Parameters
    ----------
    matrix:
        The matrix to apply householder reflections upon.
    overwrite_original:
        Whether to overwrite on the original provided matrix or not. Defaults to False.

    Returns
    -------
        Householder transformed matrix.
    """

    # taken from https://core.ac.uk/download/pdf/215673996.pdf
    matrix_ = matrix if overwrite_original else Matrix(matrix.elements[:])

    # check if there's only a single element in the matrix or not, it'll be dealt as a special case
    cond = matrix_.n_rows > 1 and matrix_.n_cols > 1

    # keeping it in if/else because it is more readable.
    if cond:
        a0_ = Matrix(matrix_.t.elements[0]).t
        a_mag = vector_mag(a0_)

        d1_ = -a_mag if a0_[0][0] > 0 else a_mag
        w11 = (a0_[0] - d1_).elements[0]

        # special case, when w11 = 0 meaning the sub-matrix only has a single element.
        if w11 == 0:
            return [matrix_, Matrix([-1])]

        v_vector = (Matrix([w11] + a0_.t.elements[1:]) / sqrt(-2 * w11 * d1_)).t

        # create a list to hold the householder transformations
        household_a = [0] * matrix.n_cols
        household_a[0] = Matrix([d1_] + [0] * (a0_.n_rows - 1)).elements

        # I - 2*((v*v.t)/(v.t*v))
        household_h = identity_matrix(matrix_.n_rows) - 2 * ((v_vector * v_vector.t) / (v_vector.t * v_vector))

        # calculate 1:n household transformations
        for i in range(1, matrix_.n_cols):
            f2 = 2 * v_vector.t * matrix_.t[i].t
            household_a[i] = (matrix_.t[i].t - f2 * v_vector).t.elements
    else:
        household_h = Matrix([1])
        household_a = (-1 * matrix_).elements[0]

    return [Matrix(household_a).t, household_h]


def populate_identity_matrix(sub_matrix: Matrix, n_rows: int, n_cols: int, s_row: int, s_col: int) -> Matrix:
    """
    Populate the identity matrix with given sub-matrices.

    Parameters
    ----------
    sub_matrix:
        The sub-matrix to populate the identity matrix into.
    n_rows:
        The number of rows in the parent matrix.
    n_cols:
        The number of columns in the parent matrix.
    s_row:
        The starting row for insertion of the sub-matrix.
    s_col:
        The starting column for insertion of the sub-matrix.

    Returns
    -------
        Identity matrix, populated with the provided sub-matrix.
    """

    # create a simple identity matrix
    populated_identity_matrix = identity_matrix(n_rows, n_cols).elements

    for i in range(sub_matrix.n_rows):
        # the sub-matrix can be populated within identity matrix easily if sub-matrix is not a single element matrix.
        if sub_matrix.n_rows > 1 and sub_matrix.n_cols > 1:
            populated_identity_matrix[s_row:][i][s_col:] = sub_matrix.elements[i]
        # if the sub-matrix only has a single element, treat it as a special case.
        else:
            populated_identity_matrix[-1][-1] = sub_matrix.elements[0]

    return Matrix(populated_identity_matrix)


def qr_decomposition_householder(matrix: Matrix, overwrite_original: bool = False) -> LMat:
    """
    Performs QR decomposition of matrix using Householder's reflections.

    Parameters
    ----------
    matrix:
        The matrix to calculate QR decomposition for.
    overwrite_original:
        Whether to overwrite on the original provided matrix or not. Defaults to False.

    Returns
    -------
        QR decomposition of the matrix.
    """

    matrix_ = matrix if overwrite_original else Matrix(matrix.elements[:])
    household_matrices: list = [0] * matrix_.n_rows
    h_matrices: list = [0] * matrix_.n_rows

    household_matrices[0], h_matrices[0] = householder_reduction(matrix_)

    if matrix_.n_rows == 2:
        return [household_matrices[0], h_matrices[0]]

    temp_ = []
    for j in range(1, matrix_.n_rows):
        temp_.append(Matrix([i[1:] for i in household_matrices[j - 1].elements[1:]]))
        # have to do separate variable insertion in lists, doesn't work with simultaneous appending
        reductions_ = householder_reduction(temp_[j - 1])
        household_matrices[j] = reductions_[0]
        h_matrices[j] = reductions_[1]

    for i in range(1, matrix_.n_rows):
        f1 = matrix_.n_rows - h_matrices[i].n_rows
        h_matrices[i] = populate_identity_matrix(h_matrices[i], matrix_.n_rows, matrix_.n_rows, f1, f1)

    q_matrix = reduce(mul, h_matrices)
    r_matrix = reduce(mul, h_matrices[::-1] + [matrix_])

    q_matrix = tolerance_to_zeros(q_matrix)
    r_matrix = tolerance_to_zeros(r_matrix)

    return [q_matrix, r_matrix]


def tolerance_to_zeros(matrix: Matrix, tolerance: IFloat = TOLERANCE) -> Matrix:
    """
    Replaces the values in the matrices below tolerance to 0s.

    Parameters
    ----------
    matrix:
        The matrix in which the values are to be replaced with zeros.
    tolerance:
        The tolerance level below which the values are reduced to zero. Default is 1e-8.

    Returns
    -------
        The matrix with 0s instead of very small values.
    """

    for i, row in enumerate(matrix.elements):
        for j, element in enumerate(row):
            if abs(element) < tolerance:
                matrix[i][j] = 0

    return matrix


def eigen_qr(matrix: Matrix, n_iter: int = 100) -> LMat:
    identity_ = identity_matrix(matrix.n_rows, matrix.n_cols)

    for _ in range(n_iter):
        temp_ = qr_decomposition_householder(matrix)
        identity_ *= temp_[0]
        matrix = mul(*temp_[::-1])

    return [identity_, matrix.diagonal()]
