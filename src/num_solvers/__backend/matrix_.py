"""Matrix_

This module provides a helpful functionality to display the matrices in fractions via,

- :class:`InFractions`

Additionally, it provides other functionalities for matrices via,

- reduce_to_zeros: Reduces the values below tolerance level to 0.
- remove_zeroed_columns: Remove the columns/rows below the all zero column/rows. Mainly for usage in QR decomposition.
- round_matrix_: Rounds the values of matrices.
- map_to_matrix: Maps a certain function to the entire matrix.
- copy_matrix: Copies or creates a deep-copy of the given matrix.

Created on Jan 10 00:01:13 2024
"""

__all__ = ['map_to_matrix', 'remove_zeroed_columns', 'round_matrix_', 'reduce_to_zeros', 'cholesky_sanity_check',
           'matrix_copy', 'upper_lower_matrices', 'gauss_methods_sanity_check', 'upper_diagonal_lower_matrices']

from copy import deepcopy

from umatrix.matrix import Matrix, identity_matrix, null_matrix

from .core_helpers_ import round_list_
from .errors_ import NonSymmetricMatrix, NotPositiveDefinite
from .. import FList, Func, IFloat, LMat, MatOrLList, N_DECIMAL, TOLERANCE


def matrix_copy(matrix: MatOrLList, overwrite: bool = False) -> Matrix:
    """
    Copy or make a deep-copy of the given matrix.

    Parameters
    ----------
    matrix:
        The matrix to make the copy of.
    overwrite:
        Whether to overwrite the original matrix or not.

    Returns
    -------
    Matrix:
        The copied matrix instance.
    """

    temp_ = Matrix(matrix) if not isinstance(matrix, Matrix) else matrix
    return Matrix(deepcopy(temp_.elements[:])) if overwrite else temp_


def reduce_to_zeros(matrix: MatOrLList, tolerance: IFloat = TOLERANCE) -> Matrix:
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

    matrix = matrix_copy(matrix)

    for i, row in enumerate(matrix.elements):
        for j, element in enumerate(row):
            if abs(element) < tolerance:
                matrix[i][j] = 0

    return matrix


def remove_zeroed_columns(matrix: MatOrLList) -> Matrix:
    """
    Removes the extra columns from the final Q and R decomposed matrices.

    Parameters
    ----------
    matrix:
        The decomposed Q or R matrix.

    Returns
    -------
        Matrix with removed columns after all 0s.
    """

    matrix = matrix_copy(matrix)

    new_, all_zeros = [0] * matrix.n_rows, []

    for i, elem in enumerate(matrix.elements):
        if not all(j == 0 for j in elem):
            new_[i] = elem
        else:
            new_[i] = elem
            all_zeros.append(i)

    if bool(all_zeros):
        if len(new_) == all_zeros[0]:
            return Matrix(new_)
        else:
            cond = len(matrix) - all_zeros[0]
            return Matrix(new_[:-cond])
    else:
        return Matrix(new_)


def round_matrix_(matrix: MatOrLList, n_decimal: int = N_DECIMAL) -> Matrix:
    """
    Maps the round function to a matrix.

    Parameters
    ----------
    matrix:
        Given matrix to round off.
    n_decimal:
        Number of decimal places to round off to.

    Returns
    -------
        Rounded off matrix.
    """

    matrix = matrix_copy(matrix)

    return Matrix(round_list_(matrix.elements, n_decimal))


def map_to_matrix(matrix: MatOrLList, function: Func):
    """
    Apply a given function element-wise to a matrix.

    Parameters
    ----------
    matrix:
        The matrix to be mapped.
    function:
        A function that takes a single float as input and returns a float.

    Returns
    -------
    Matrix
        A new matrix where the function has been applied element-wise.
    """

    matrix = matrix_copy(matrix)

    for i in range(matrix.n_rows):
        for j in range(matrix.n_cols):
            matrix[j][j] = function(matrix[i][j])

    return matrix


def cholesky_sanity_check(matrix: MatOrLList) -> Matrix:
    """
    Performs sanity check for Cholesky decomposition.

    Parameters
    ----------
    matrix:
        The matrix to perform cholesky decomposition on.

    Raises
    ------
    NonSymmetricMatrix:
        If the matrix is not symmetric.
    NotPositiveDefinite:
        If the matrix is not positive definite.
    """

    matrix = matrix_copy(matrix)

    if not matrix.is_symmetric():
        raise NonSymmetricMatrix('The matrix is not symmetric. Can not perform Cholesky decomposition.')

    if not matrix.is_positive_definite():
        raise NotPositiveDefinite('The matrix is not positive definite. Can not perform Cholesky decomposition.')

    return matrix


def upper_diagonal_lower_matrices(matrix: Matrix) -> LMat:
    """
    Decomposes a matrix into its upper, diagonal and lower triangular parts.

    Parameters
    ----------
    matrix:
        The input matrix.

    Returns
    -------
    LMat
        Lower, Diagonal and Upper triangular matrix.
    """

    n_rows, n_cols, matrix_l = matrix.n_rows, matrix.n_cols, null_matrix(matrix.n_rows)

    for i in range(n_rows):
        for j in range(n_cols):
            if i > j:
                matrix_l[i][j] = matrix[i][j]

    matrix_d = matrix.diagonal_of_matrix()

    return [matrix_l, matrix_d, matrix - matrix_l - matrix_d]


def upper_lower_matrices(matrix: Matrix) -> LMat:
    """
    Decomposes a matrix into its upper and lower triangular parts.

    Parameters
    ----------
    matrix:
        The input matrix.

    Returns
    -------
    LMat
        Lower and Upper triangular matrix.
    """

    n_rows, n_cols, matrix_l = matrix.n_rows, matrix.n_cols, identity_matrix(matrix.n_rows)

    for i in range(n_rows):
        for j in range(n_cols):
            if i >= j:
                matrix_l[i][j] = matrix[i][j]

    return [matrix_l, matrix - matrix_l]


def gauss_methods_sanity_check(matrix: Matrix, solution: FList, initial_guess: FList):
    def is_diagonally_dominant():
        matrix_diagonal = matrix_.diagonal()
        temp_ = [matrix_diagonal[j] >= (sum(matrix_[j]) - matrix_diagonal[j]) for j in range(matrix_.n_rows)]

        return all(x for x in temp_)

    matrix_, solution, initial_guess = matrix_copy(matrix), Matrix(solution).t, Matrix(initial_guess).t

    if not matrix_.is_symmetric():
        print('The matrix is not symmetric.')

    if not matrix_.is_positive_definite():
        print('The matrix is not positive definite.')

    if not is_diagonally_dominant():
        print('The system of equations is not diagonally dominant. The solutions might not be correct.')

    return matrix_, solution, initial_guess
