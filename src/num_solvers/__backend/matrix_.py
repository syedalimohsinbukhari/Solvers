"""Matrix_

This module provides a helpful functionality to display the matrices in fractions via,

- :class:`InFractions`

Additionally, it provides other functionalities for matrices via,

- reduce_to_zeros: Reduces the values below tolerance level to 0.
- remove_zeroed_columns: Remove the columns/rows below the all zero column/rows. Mainly for usage in QR decomposition.
- round_matrix_: Rounds the values of matrices.

Created on Jan 10 00:01:13 2024
"""

__all__ = ['remove_zeroed_columns', 'round_matrix_', 'reduce_to_zeros', 'cholesky_sanity_check',
           'upper_lower_matrices', 'gauss_methods_sanity_check', 'upper_diagonal_lower_matrices',
           'populate_identity_matrix']

from umatrix.matrix import Matrix, identity_matrix, matrix_copy, null_matrix

from .core_helpers_ import round_list_
from .errors_ import NonSymmetricMatrix, NotPositiveDefinite
from .. import FList, IFloat, LMat, MatOrLList, N_DECIMAL, TOLERANCE


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


def strassen_mat_mul(matrix_a: Matrix, matrix_b: Matrix):
    m1 = (matrix_a[0][0] + matrix_a[1][1]) * (matrix_b[0][0] + matrix_b[1][1])
    m2 = (matrix_a[1][0] + matrix_a[1][1]) * matrix_b[0][0]
    m3 = matrix_a[0][0] * (matrix_b[0][1] - matrix_b[1][1])
    m4 = matrix_a[1][1] * (matrix_b[1][0] - matrix_b[0][0])
    m5 = (matrix_a[0][0] + matrix_a[0][1]) * matrix_b[1][1]
    m6 = (matrix_a[1][0] - matrix_a[0][0]) * (matrix_b[0][0] + matrix_b[0][1])
    m7 = (matrix_a[0][1] - matrix_a[1][1]) * (matrix_b[1][0] + matrix_b[1][1])

    c11 = m1 + m4 - m5 + m7
    c12 = m3 + m5
    c21 = m2 + m4
    c22 = m1 - m2 + m3 + m6

    return Matrix([[c11, c12], [c21, c22]])
