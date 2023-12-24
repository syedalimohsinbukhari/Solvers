"""Created on Dec 20 23:21:44 2023"""

from cmath import sqrt
from typing import List

import numpy as np

TOLERANCE = 1e-8


def generate_random_polynomial(degree: int, low: int = -10, high: int = 10) -> List[float]:
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


def laguerre_method(polynomial: List, x_0: float, degree_of_polynomial: int, tolerance: float = TOLERANCE,
                    get_full_result: bool = False):
    p_val, p_der = np.polyval, np.polyder
    root_ = [x_0]

    while True:
        poly = p_val(polynomial, x_0)

        if abs(poly) < tolerance:
            break

        poly_derivative = p_val(p_der(polynomial), x_0)
        poly_double_derivative = p_val(p_der(polynomial, 2), x_0)

        g_ = poly_derivative / poly
        h_ = g_**2 - (poly_double_derivative / poly)

        sqrt_in_denominator = sqrt((degree_of_polynomial - 1) * (degree_of_polynomial * h_ - g_**2))
        p1 = max([g_ + sqrt_in_denominator, g_ - sqrt_in_denominator], key=abs)

        x_0 = x_0 - (degree_of_polynomial / p1)

        root_.append(x_0)

    return root_ if get_full_result else root_[-1]


def segmented_roots(polynomial: List, x_start: float, x_end: float, degree_of_polynomial: int, num_segments: int = 100,
                    tolerance: float = TOLERANCE):
    x_values = np.linspace(x_start, x_end, num_segments)
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

    return new_roots
