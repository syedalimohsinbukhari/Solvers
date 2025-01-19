"""
Singular Value Decomposition using Power Iteration

This module provides an implementation of the Singular Value Decomposition (SVD) using the Power Iteration method.
The SVD decomposes a matrix into the product of three matrices: U, Sigma, and V^T.

References:
    https://www.cs.yale.edu/homes/el327/datamining2013aFiles/07_singular_value_decomposition.pdf

Created on Jan 18 00:48:38 2024
"""

import random

from umatrix.matrix import Matrix, identity_matrix, vector_mag

from .. import LMat, TOLERANCE


# works for mxn where m=n
# works for mxn where m=n-1

def svd_power_method(matrix: Matrix) -> LMat:
    """
    Compute the Singular Value Decomposition (SVD) of a given matrix using Power Iteration.

    Parameters
    ----------
    matrix : Matrix
        The input matrix for which the SVD is to be computed.

    Returns
    -------
    List[Matrix]
        A list containing three matrices: U, Sigma, and V^T, such that matrix = U * Sigma * V^T.

    Notes
    -----
    This implementation uses the Power Iteration method to compute the singular value decomposition.
    It works for matrices with dimensions m x n where m=n or m=n-1.

    References
    ----------
    https://www.cs.yale.edu/homes/el327/datamining2013aFiles/07_singular_value_decomposition.pdf

    """

    def uv_difference(u_or_v_list, difference):
        temp_ = []
        for element in u_or_v_list[:matrix.n_cols - difference]:
            temp_.append(element.t.elements[:matrix.n_cols - difference])

        return temp_

    def sigma_difference(sigma_list, difference):
        temp_ = identity_matrix(len(sigma))
        for element in range(temp_.n_cols):
            temp_[element][element] = sigma_list[element]

        temp_ = Matrix([element[:matrix.n_cols - difference] for element in temp_.elements])

        return temp_

    # get the original n_rows and n_cols
    nr, nc = matrix.n_rows, matrix.n_cols

    # if original nr > nc, transpose the matrix, it's easier to deal with nc > nr
    if nr > nc:
        matrix = matrix.t

    # make a copy
    matrix_ = Matrix(matrix.t.elements[:])

    # add tolerance to excessive columns
    if matrix_.n_rows > matrix_.n_cols:
        for v in matrix_.elements:
            v.extend([TOLERANCE])

    # initialize u, sigma, and u
    v: list = [0] * matrix_.n_rows
    sigma: list = [0] * matrix_.n_rows
    u: list = [0] * matrix_.n_cols

    # iterate through the number of rows using power-iteration method
    for j in range(matrix_.n_rows):
        x0 = Matrix([random.random() for _ in range(matrix_.n_cols)]).t
        x: list = [0] * 10
        x[0] = x0

        for i in range(1, 10):
            x[i] = matrix_.t * matrix_ * x[i - 1]

        v[j] = x[-1] / vector_mag(x[-1])
        sigma[j] = vector_mag(matrix_ * v[j])
        u[j] = (matrix_ * v[j]) / sigma[j]

        matrix_ -= u[j] * sigma[j] * v[j].t

    # re-calibrate the u, v, and sigma matrices.
    u_matrix = uv_difference(u, 0)

    if nr != nc:
        # re-calibrate the u, v, and sigma matrices.
        v_matrix = uv_difference(v, 1)
        sigma_matrix = sigma_difference(sigma, 1)
    else:
        v_matrix = uv_difference(v, 0)
        sigma_matrix = sigma_difference(sigma, 0)

    # change list of list to matrices
    u_matrix = Matrix(u_matrix).t
    v_matrix = Matrix(v_matrix)
    sigma_matrix = sigma_matrix

    # return the result, u, sigma, and v.t
    # the user only needs to multiply the given matrices to get the original one.
    if nc > nr:
        v_matrix, u_matrix = u_matrix, v_matrix
        return [u_matrix.t, sigma_matrix.t, v_matrix.t]
    else:
        return [u_matrix, sigma_matrix, v_matrix]
