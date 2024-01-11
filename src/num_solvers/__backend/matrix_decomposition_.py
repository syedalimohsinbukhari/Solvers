"""Created on Jan 12 00:59:01 2024"""

from .extra_ import round_list_
from .. import LLList, LList, N_DECIMAL


def __lower_upper_matrices(len_mat: int) -> LLList:
    """
    Generate lower and upper matrices with zeros.

    Parameters
    ----------
    len_mat:
        Size of the matrices.

    Returns
    -------
        List containing the lower and upper matrices filled with zeros.
    """

    lower_matrix = [[0] * len_mat for _ in range(len_mat)]
    upper_matrix = [[0] * len_mat for _ in range(len_mat)]

    return [lower_matrix, upper_matrix]


def __round_matrices(lower_matrix: LLList, upper_matrix: LLList, n_decimal: int = N_DECIMAL) -> LLList:
    """
    Round all values in the lower and upper matrices.

    Parameters
    ----------
    lower_matrix:
        Lower matrix to be rounded.
    upper_matrix:
        Upper matrix to be rounded.
    n_decimal:
        Number of decimals to round the matrices to (default is N_DECIMAL).

    Returns
    -------
        List containing the rounded lower and upper matrices.
    """

    lower_matrix = [round_list_(row, n_decimal) for row in lower_matrix]
    upper_matrix = [round_list_(row, n_decimal) for row in upper_matrix]

    return [lower_matrix, upper_matrix]


def __transpose(matrix: LList) -> LList:
    """
    Transpose a given matrix.

    Parameters
    ----------
    matrix:
        Matrix to be transposed.

    Returns
    -------
        Transposed matrix.
    """

    n_dim = len(matrix)
    return [[matrix[row][col] for row in range(n_dim)] for col in range(n_dim)]


def __is_symmetric(matrix: LList) -> bool:
    """Checks if the given matrix is symmetric or not.

    Parameters
    ----------
    matrix:
        The matrix to check.

    Returns
    -------
        True if the matrix is symmetric, else False.
    """

    return matrix == __transpose(matrix)


def __is_positive_definite(matrix: LList) -> bool:
    """Checks if the matrix is positive definite or not.

    Parameters
    ----------
    matrix:
        The matrix to check.

    Returns
    -------
        True if the matrix is positive definite, else False.

    """

    n_dimensions = len(matrix)

    if len(matrix) != len(matrix[0]):
        return False

    if not __is_symmetric(matrix):
        return False

    for i in range(1, n_dimensions + 1):
        minor = [row[:i] for row in matrix[:i]]
        determinant = sum((minor[j][j] for j in range(i)))
        if determinant <= 0:
            return False

    return True
