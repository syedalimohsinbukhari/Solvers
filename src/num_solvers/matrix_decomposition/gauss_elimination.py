"""Created on Jan 21 00:12:13 2024"""

from .matrix import Matrix
from .. import IFloat, TOLERANCE


def gauss_elimination(matrix: Matrix, tolerance: IFloat = TOLERANCE, modify_original: bool = False) -> Matrix:
    matrix_ = matrix if modify_original else Matrix(matrix.elements[:])
    num_rows, num_cols = matrix_.n_rows, matrix.n_cols

    for i in range(min(num_rows, num_cols)):
        d_element = matrix_[i][i]

        cond1 = abs(d_element) > tolerance

        if cond1:
            matrix_[i] = [elem / d_element
                          if cond1 else 0.0
                          for elem in matrix_[i]]

        for j in range(num_rows):
            if i != j:
                factor = matrix_[j][i]
                matrix_[j] = [a - factor * b
                              if abs(a - factor * b) > tolerance else 0.0
                              for a, b in zip(matrix_[j], matrix_[i])]

    return matrix_