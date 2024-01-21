"""Created on Jan 21 15:05:20 2024"""

import functools
import operator
from _operator import mul
from collections import OrderedDict
from math import sqrt

from src.num_solvers import TOLERANCE
from src.num_solvers.matrix_decomposition import LMat, Matrix
from src.num_solvers.matrix_decomposition.matrix import identity_matrix, vector_mag


def householder_reduction(matrix: Matrix, overwrite_original: bool = False) -> LMat:
    # taken from https://core.ac.uk/download/pdf/215673996.pdf
    matrix_ = matrix if overwrite_original else Matrix(matrix.elements[:])

    cond = matrix_.n_rows > 1 and matrix_.n_cols > 1

    if cond:
        a0 = Matrix(matrix_.t.elements[0]).t
        a_mag = vector_mag(a0)
    else:
        a0 = Matrix([[matrix_.t.elements]]).t
        a_mag = a0.elements[0][0]

    d1 = -a_mag if a0[0][0] > 0 else a_mag
    w11 = (a0[0] - d1).elements[0]

    if w11 == 0:
        return [matrix_, Matrix([-1])]

    f1 = sqrt(-2 * w11 * d1)
    v = (Matrix([w11] + a0.t.elements[1:]) / f1).t

    household_a = [0] * matrix.n_cols
    household_a[0] = Matrix([d1] + [0] * (a0.n_rows - 1)).elements

    household_h = identity_matrix(matrix_.n_rows) - 2 * ((v * v.t) / (v.t * v))

    for i in range(1, matrix_.n_cols):
        f2 = 2 * v.t * matrix_.t[i].t
        household_a[i] = (matrix_.t[i].t - f2 * v).t.elements

    return [Matrix(household_a).t, household_h]


def populate_identity_matrix(matrix: Matrix, n_rows: int, n_cols: int, s_row: int, s_col: int):
    populated_identity_matrix = identity_matrix(n_rows, n_cols).elements
    for i in range(matrix.n_rows):
        if matrix.n_rows > 1 and matrix.n_cols > 1:
            populated_identity_matrix[s_row:][i][s_col:] = matrix.elements[i]
        else:
            populated_identity_matrix[-1][-1] = matrix.elements[0]

    return Matrix(populated_identity_matrix)


def qr_decomposition_householder(matrix: Matrix, overwrite_original: bool = False):
    matrix_ = matrix if overwrite_original else Matrix(matrix.elements[:])
    household_matrices = [0] * matrix_.n_rows
    h_matrices = [0] * matrix_.n_rows

    household_matrices[0], h_matrices[0] = householder_reduction(matrix_)

    temp_ = []
    for j in range(1, matrix_.n_rows):
        temp_.append(Matrix([i[1:] for i in household_matrices[j - 1].elements[1:]]))
        f1 = householder_reduction(temp_[j - 1])
        household_matrices[j] = f1[0]
        h_matrices[j] = f1[1]

    for i in range(1, matrix_.n_rows):
        f1 = matrix_.n_rows - h_matrices[i].n_rows
        h_matrices[i] = populate_identity_matrix(h_matrices[i],
                                                 matrix_.n_rows,
                                                 matrix.n_cols,
                                                 f1,
                                                 f1)

    q = functools.reduce(operator.mul, h_matrices)
    r = functools.reduce(operator.mul, h_matrices[::-1] + [matrix_])

    return OrderedDict({'Q': get_zeros(q), 'R': get_zeros(r)})


def get_zeros(matrix: Matrix) -> Matrix:
    for i, row in enumerate(matrix.elements):
        for j, element in enumerate(row):
            if abs(element) < TOLERANCE:
                matrix[i][j] = 0

    return matrix


def eeigen_qr(matrix: Matrix, n_iter: int = 100) -> LMat:
    identity_ = identity_matrix(matrix.n_rows, matrix.n_cols)

    for _ in range(n_iter):
        temp_ = qr_decomposition_householder(matrix)
        identity_ *= temp_[0]
        matrix = mul(*temp_[::-1])

    return [identity_, matrix.diagonal()]
