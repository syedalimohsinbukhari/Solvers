"""Polynomial root finding algorithms

This module provides functionality to find roots of given polynomials. The roots are determined via,

- laguerre_method

In addition, the user can specify a range in which all roots are to be determined for the given polynomial via,

- segmented_roots

This module also provides a handy function to determine random polynomials, via

- generate_random_polynomial

Created on Dec 20 23:21:44 2023
"""

from cmath import sqrt
from typing import List

import numpy as np

from .. import FList, IFloat, LList, TOLERANCE


# TODO: See if numpy can be removed


def laguerre_method(polynomial: FList, x_guess: IFloat, degree_of_polynomial: int, tolerance: IFloat = TOLERANCE,
                    get_full_result: bool = False):
    """
    Use Laguerre method to find the roots of the given polynomial.

    Parameters
    ----------
    polynomial:
        The polynomial as list of coefficients.
    x_guess:
        Initial guess for the root.
    degree_of_polynomial:
        Degree of the polynomial.
    tolerance:
        Tolerance for the result. Default is TOLERANCE.
    get_full_result:
        If True, gives all the guesses. If false, only gives the root of polynomial.

    Returns
    -------
        Root of the given polynomial.
    """

    p_val, p_der = np.polyval, np.polyder
    root_ = [x_guess]

    while True:
        poly = p_val(polynomial, x_guess)

        if abs(poly) < tolerance:
            break

        poly_derivative = p_val(p_der(polynomial), x_guess)
        poly_double_derivative = p_val(p_der(polynomial, 2), x_guess)

        g_ = poly_derivative / poly
        h_ = g_**2 - (poly_double_derivative / poly)

        sqrt_in_denominator = sqrt((degree_of_polynomial - 1) * (degree_of_polynomial * h_ - g_**2))
        p1 = max([g_ + sqrt_in_denominator, g_ - sqrt_in_denominator], key=abs)

        x_guess = x_guess - (degree_of_polynomial / p1)

        root_.append(x_guess)

    return root_ if get_full_result else root_[-1]


def segmented_roots(polynomial: List, x_range: FList, degree_of_polynomial: int,
                    num_segments: int = 100, tolerance: IFloat = TOLERANCE) -> LList:
    """
    Segments the given interval to find all possible roots within that interval for a given polynomial.

    Parameters
    ----------
    polynomial:
        The polynomial as list of coefficients.
    x_range:
        The range in which the root of polynomial is to be determined.
    degree_of_polynomial:
        Degree of the polynomial.
    num_segments:
        Number of segments to divide the ``x_range`` into. Default is 100.
    tolerance:
        Tolerance for the result. Default is TOLERANCE.

    Returns
    -------
        List of all possible roots within the given interval.
    """

    x_values = np.linspace(x_range[0], x_range[1], num_segments)
    all_roots = []

    for x_0 in x_values:
        root = laguerre_method(polynomial, x_0, degree_of_polynomial, tolerance)

        if all(abs(root - existing_root) > tolerance for existing_root in all_roots):
            all_roots.append(root)

    unique_roots = np.fromiter(all_roots, dtype=complex)

    new_roots = []
    for i, root_ in enumerate(unique_roots):
        if np.allclose(np.imag(root_), 0):
            new_roots.append(np.real(root_))
        else:
            new_roots.extend([root_, np.conj(root_)])

    new_roots = list(map(lambda x: round(x, 8), new_roots))

    return new_roots


def generate_random_polynomial(degree: int, low: int = -10, high: int = 10) -> FList:
    """
    Generate random coefficients for a polynomial.

    Parameters
    ----------
    degree : int
        The degree of the polynomial.
    low : int, optional
        The lower bound for the random coefficients. Default is -10.
    high : int, optional
        The upper bound for the random coefficients. Default is 10.

    Returns
    -------
    List[float]
        A list of coefficients representing the random polynomial.
    """

    coefficients = np.random.uniform(low, high, degree + 1)

    while coefficients[0] == 0:
        coefficients[0] = np.random.uniform(low, high)

    return list(coefficients)
