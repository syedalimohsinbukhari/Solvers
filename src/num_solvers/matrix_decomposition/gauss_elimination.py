"""Gauss elimination

This module provides functionality to apply gauss elimination process on given matrices, via,

- gauss_elimination


Created on Jan 21 00:12:13 2024
"""

__all__ = ['gauss_elimination']

from umatrix.matrix import Matrix, matrix_copy

from .. import IFloat, MatOrLList, TOLERANCE


def gauss_elimination(matrix: MatOrLList, tolerance: IFloat = TOLERANCE, modify_original: bool = False) -> Matrix:
    """
    Performs Gaussian elimination on the given matrix.

    Parameters
    ----------
    matrix:
        The matrix to perform Gaussian elimination on.
    tolerance:
        Tolerance for convergence. Default is TOLERANCE.
    modify_original:
        Whether to modify the original given matrix or not.

    Returns
    -------
    Matrix:
        The modified matrix with Gaussian elimination applied.
    """

    matrix_ = matrix_copy(matrix, modify_original)

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
