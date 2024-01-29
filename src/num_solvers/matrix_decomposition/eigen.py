"""Created on Jan 28 15:42:03 2024"""

from math import sqrt

from umatrix.matrix import Matrix, identity_matrix

from .. import FList


def eigenvalues_2x2(matrix: Matrix) -> FList:
    """
    Compute the eigenvalues of a 2x2 matrix.

    Parameters
    ----------
    matrix : np.ndarray
        The input 2x2 matrix.

    Returns
    -------
    List[float]
        A list containing the two eigenvalues.
    """

    trace = matrix.trace
    determinant = matrix.determinant()

    num1 = trace + sqrt(trace**2 - 4 * determinant)
    num2 = trace - sqrt(trace**2 - 4 * determinant)

    return [num1 / 2, num2 / 2]


def eigenvectors_2x2(matrix: Matrix):
    """
    Compute the eigenvectors of a 2x2 matrix.

    Parameters
    ----------
    matrix : np.ndarray
        The input 2x2 matrix.

    Returns
    -------
    Tuple[List[float], List[np.ndarray]]
        A tuple containing a list of eigenvalues
        and a list of corresponding eigenvectors.
    """

    eigenvalues = eigenvalues_2x2(matrix)
    eigenvectors = [matrix - identity_matrix(2, value=value) for value in eigenvalues[::-1]]

    scaled_eigenvectors = [(eigenvector.t[0] / eigenvector.t[0][-1]).elements for eigenvector in eigenvectors]

    return eigenvalues, Matrix(scaled_eigenvectors)
