"""Householder transformations

This module provides functionality to apply Householder transformation to a given matrix, via,

- householder_reduction: Applies the Householder reduction to the columns of the given matrix.
- qr_decomposition_householder: Decomposes the given matrix into Q and R matrices via householder reductions.

Created on Jan 21 15:05:20 2024
"""

__all__ = ['householder_reduction', 'qr_decomposition_householder']

from functools import reduce
from math import sqrt
from operator import mul

from umatrix.matrix import Matrix, identity_matrix, vector_mag

from .. import LMat, MatOrLList
from ..__backend.matrix_ import populate_identity_matrix, reduce_to_zeros


# taken from https://core.ac.uk/download/pdf/215673996.pdf
def householder_reduction(matrix: MatOrLList, overwrite_original: bool = False) -> LMat:
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

    def calculate_householder_parameters():
        a_mag = vector_mag(a0_)

        d1 = -a_mag if a0_[0][0] > 0 else a_mag
        w11 = (a0_[0] - d1).elements[0]

        # special case, when w11 = 0 meaning the sub-matrix only has a single element.
        if w11 == 0:
            return [matrix_, Matrix([-1])]

        v_vector = (Matrix([w11] + a0_.t.elements[1:]) / sqrt(-2 * w11 * d1)).t

        return d1, v_vector

    if not isinstance(matrix, Matrix):
        matrix = Matrix(matrix)

    matrix_: Matrix = matrix if overwrite_original else Matrix(matrix.elements[:])

    # check if there's only a single element in the matrix or not, it'll be dealt as a special case
    cond = matrix_.n_rows > 1 and matrix_.n_cols > 1

    # keeping it in if/else because it is more readable.
    if cond:
        a0_ = Matrix(matrix_.t.elements[0]).t

        d1_, v_vector_ = calculate_householder_parameters()

        # create a list to hold the householder transformations
        household_a = [0] * matrix.n_cols
        household_a[0] = Matrix([d1_] + [0] * (a0_.n_rows - 1)).elements

        # I - 2*((v*v.t)/(v.t*v))
        household_h = identity_matrix(matrix_.n_rows) - 2 * ((v_vector_ * v_vector_.t) / (v_vector_.t * v_vector_))

        # calculate 1:n household transformations
        for i in range(1, matrix_.n_cols):
            f1 = 2 * v_vector_.t * matrix_.t[i].t
            household_a[i] = (matrix_.t[i].t - f1 * v_vector_).t.elements
    else:
        if matrix_.n_rows == 1:
            household_h = Matrix([1])
            household_a = (-1 * matrix_).elements[0]
        else:
            a0_ = matrix_
            d1_, v_vector_ = calculate_householder_parameters()

            household_h = identity_matrix(matrix_.n_rows) - 2 * ((v_vector_ * v_vector_.t) / (v_vector_.t * v_vector_))
            household_a = (household_h * a0_).elements

    return [Matrix(household_a).t, household_h]


def qr_decomposition_householder(matrix: MatOrLList, overwrite_original: bool = False) -> LMat:
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

    if not isinstance(matrix, Matrix):
        matrix = Matrix(matrix)

    matrix_: Matrix = matrix if overwrite_original else Matrix(matrix.elements[:])

    if matrix_.n_rows > matrix_.n_cols:
        remove_ = True
        cond = matrix_.n_cols
    else:
        remove_ = False
        cond = matrix_.n_rows

    household_matrices: list = [0] * cond
    h_matrices: list = [0] * cond

    household_matrices[0], h_matrices[0] = householder_reduction(matrix_)

    if matrix_.n_rows == 2:
        return [household_matrices[0], h_matrices[0]]

    temp_ = []
    for j in range(1, cond):
        temp_.append(Matrix([i[1:] for i in household_matrices[j - 1].elements[1:]]))
        # have to do separate variable insertion in lists, doesn't work with simultaneous appending
        reductions_ = householder_reduction(temp_[j - 1])
        household_matrices[j] = reductions_[0]
        h_matrices[j] = reductions_[1]

    for i in range(1, cond):
        f1 = matrix_.n_rows - h_matrices[i].n_rows
        h_matrices[i] = populate_identity_matrix(h_matrices[i], matrix_.n_rows, matrix_.n_rows, f1, f1)

    q_matrix = reduce(mul, h_matrices)
    r_matrix = reduce(mul, h_matrices[::-1] + [matrix_])

    q_matrix = reduce_to_zeros(q_matrix)
    r_matrix = reduce_to_zeros(r_matrix)

    if remove_:
        q_matrix = Matrix(q_matrix.t.elements[:-1]).t
        r_matrix = Matrix(r_matrix.elements[:-1])

    return [q_matrix, r_matrix]

# def eigen_qr(matrix: MatOrLList, n_iter: int = 100) -> LMat:
#     if not isinstance(matrix, Matrix):
#         matrix = Matrix(matrix)
#
#     identity_ = identity_matrix(matrix.n_rows, matrix.n_cols)
#
#     for _ in range(n_iter):
#         temp_ = qr_decomposition_householder(matrix)
#         identity_ *= temp_[0]
#         matrix = mul(*temp_[::-1])
#
#     return [identity_, matrix.diagonal()]
