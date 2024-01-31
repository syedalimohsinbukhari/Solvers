"""Extra

This module provides some extra functionalities related to matrix operations, via,

- characteristic_polynomial: Provides characteristic polynomial for the given matrix.
- faddeev_le_verrier: The parent function for characteristic polynomial calculation of the matrix.
- eigenvalues_2x2: Provides eigen values for 2x2 matrices.
- eigenvectors_2x2: Provides eigen vectors for 2x2 matrices.

Created on Jan 31 12:45:37 2024
"""

__all__ = ['characteristic_polynomial', 'eigenvectors_2x2', 'eigenvalues_2x2']

from math import sqrt

from custom_inherit import doc_inherit
from umatrix.matrix import Matrix, identity_matrix

from .. import DOC_STYLE, FList, LList, MatOrLList


def faddeev_le_verrier(matrix: MatOrLList) -> FList:
    """
    Calculates the characteristic polynomial of the given matrix using Faddeev Le'Verrier method.

    Parameters
    ----------
    matrix:
        Given matrix to calculate characteristic polynomial for.

    Returns
    -------
    FList:
        List of polynomial coefficients corresponding to the characteristic polynomial for the given matrix.
    """

    new_matrix_: list = [0] * (matrix.n_cols + 1)
    coefficient: list = [0] * (matrix.n_cols + 1)

    new_matrix_[0], coefficient[0] = identity_matrix(matrix.n_rows), 1

    for i in range(1, matrix.n_rows + 1):
        temp_ = matrix * new_matrix_[i - 1]
        coefficient[i] = (-1 / i) * sum(temp_.diagonal().elements)
        new_matrix_[i] = temp_ + identity_matrix(matrix.n_rows, value=coefficient[i])

    return coefficient


@doc_inherit(faddeev_le_verrier, style=DOC_STYLE)
def characteristic_polynomial(matrix: MatOrLList) -> FList:
    return faddeev_le_verrier(matrix)


def eigenvalues_2x2(matrix: MatOrLList) -> FList:
    """
    Compute the eigenvalues of a 2x2 matrix.

    Parameters
    ----------
    matrix : np.ndarray
        The input 2x2 matrix.

    Returns
    -------
    FList
        A list containing the two eigenvalues.
    """

    trace = matrix.trace
    determinant = matrix.determinant()

    num1 = trace + sqrt(trace**2 - 4 * determinant)
    num2 = trace - sqrt(trace**2 - 4 * determinant)

    return [num1 / 2, num2 / 2]


def eigenvectors_2x2(matrix: MatOrLList) -> LList:
    """
    Compute the eigenvectors of a 2x2 matrix.

    Parameters
    ----------
    matrix : np.ndarray
        The input 2x2 matrix.

    Returns
    -------
    LList
        A tuple containing a list of eigenvalues and a list of corresponding eigenvectors.
    """

    eigenvalues = eigenvalues_2x2(matrix)
    eigenvectors = [matrix - identity_matrix(2, value=value) for value in eigenvalues[::-1]]

    scaled_eigenvectors = [(eigenvector.t[0] / eigenvector.t[0][-1]).elements for eigenvector in eigenvectors]

    return [eigenvalues, Matrix(scaled_eigenvectors)]
