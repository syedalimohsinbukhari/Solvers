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
from operator import sub
from typing import List

import numpy as np

from .. import FList, IFloat, LList, OptList, TOLERANCE
from ..__backend.core_helpers_ import filter_similar_values, round_list_
from ..__backend.errors_ import DegreeOfPolynomialNotCorrect


# TODO: See if numpy can be removed


def laguerre_method(polynomial: FList, degree_of_polynomial: int = -1, x_guess: IFloat = 0,
                    get_full_result: bool = False, tolerance: IFloat = TOLERANCE):
    """
    Use Laguerre method to find the roots of the given polynomial.

    Parameters
    ----------
    polynomial:
        The polynomial as list of coefficients.
    degree_of_polynomial:
        Degree of the polynomial. Defaults to -1. If polynomial's degree is -1, it is calculated inside the function.
    x_guess:
        Initial guess for the root. Defaults to 0.
    get_full_result:
        If True, gives all the guesses. If false, only gives the root of polynomial.
    tolerance:
        Tolerance for the result. Default is TOLERANCE.

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

        if len(root_) > 20 and sub(*root_[-2:][::-1]) < tolerance:
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


def segmented_roots(polynomial: List, x_guess: OptList = None,
                    num_segments: int = 100, n_decimal: int = 8, tolerance: IFloat = TOLERANCE) -> LList:
    """
    Segments the given interval to find all possible roots within that interval for a given polynomial.

    Parameters
    ----------
    polynomial:
        The polynomial as list of coefficients.
    x_guess:
        The range in which the root of polynomial is to be determined. Default is [-100, 100].
    num_segments:
        Number of segments to divide the ``x_range`` into. Default is 500.
    n_decimal:
        Number of digits to round off to. Default is 8
    tolerance:
        Tolerance for the result. Default is TOLERANCE.

    Returns
    -------
        List of all possible roots within the given interval.
    """

    x_guess = x_guess if x_guess else [-10, 50]
    degree_of_polynomial = len(polynomial) - 1

    if degree_of_polynomial:
        if len(polynomial) - 1 != degree_of_polynomial:
            raise DegreeOfPolynomialNotCorrect('The provided polynomial and degree of polynomial do not match')

    x_values = np.linspace(x_guess[0], x_guess[1], num_segments)
    all_roots = []

    for x_0 in x_values:
        root = laguerre_method(polynomial,
                               degree_of_polynomial,
                               x_0,
                               tolerance=tolerance)

        if all(abs(root - existing_root) > tolerance for existing_root in all_roots):
            all_roots.append(root)

    unique_roots = np.fromiter(all_roots, dtype=complex)

    new_roots = []
    for i, root_ in enumerate(unique_roots):
        if np.allclose(np.imag(root_), 0):
            new_roots.append(np.real(root_))
        else:
            new_roots.extend([root_, np.conj(root_)])

    new_roots = round_list_(new_roots, n_decimal)

    while len(new_roots) != degree_of_polynomial:
        new_roots = round_list_(new_roots, n_decimal - 1)
        n_decimal -= 1
        if n_decimal == degree_of_polynomial:
            break

    new_roots = filter_similar_values(new_roots)

    if len(new_roots) < degree_of_polynomial:
        print('The roots of polynomial contains repeating roots.')

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
