"""Created on Dec 24 12:56:39 2023"""

from math import sqrt
from typing import Callable, List, Union

Num = Union[float, complex]
Func = Callable[[Num], Num]
TOLERANCE = 1e-8

FLOAT_LIST = Union[List, float]


def bisection_method(function, x_0, x_n, tolerance=TOLERANCE, get_all_values=False):
    f, root = function, []

    def is_sane():
        return f(x_0) < 0 < f(x_n)

    while is_sane():
        x_new = (x_0 + x_n) / 2
        root.append(x_new)
        f_x_new = f(x_new)

        if abs(f_x_new) < tolerance:
            break

        if f_x_new < 0:
            x_0 = x_new
        else:
            x_n = x_new

    return root if get_all_values else root[-1]


def false_position_method(function, x_0, x_n, tolerance=TOLERANCE, get_all_values=False):
    f, root = function, []

    def is_sane():
        return f(x_0) < 0 < f(x_n)

    while is_sane():
        f_x0, f_x1 = f(x_0), f(x_n)

        x_new = (x_0 * f_x1 - x_n * f_x0) / (f_x1 - f_x0)
        root.append(x_new)

        f_x_new = f(x_new)

        if abs(f_x_new) < tolerance:
            break

        if f_x_new == 0:
            return root if get_all_values else root[-1]

        if f_x0 * f_x_new < 0:
            x_n = x_new
        else:
            x_0 = x_new

    return root if get_all_values else root[-1]


def regula_falsi_method(function, x_0, x_n, tolerance=TOLERANCE, get_all_values=False):
    return false_position_method(function, x_0, x_n, tolerance, get_all_values)


def secant_method(function, x_0, x_n, tolerance=TOLERANCE, get_all_values=False):
    f, root = function, []

    while True:
        x_new = x_n - f(x_n) * (x_n - x_0) / (f(x_n) - f(x_0))

        root.append(x_new)
        x_0, x_n = x_n, x_new

        f_x_new = f(x_new)
        if abs(f_x_new) < tolerance:
            break

        if f_x_new == 0:
            return root if get_all_values else root[-1]

    return root if get_all_values else root[-1]


# Taken form https://en.wikipedia.org/wiki/Muller%27s_method#Computational_example
# minor tweaking applied on variable namings for consistency
def div_diff(f_: Func, xs_: list[Num]):
    """Calculate the divided difference f[x0, x1, ...]."""
    if len(xs_) == 2:
        a, b = xs_
        return (f_(a) - f_(b)) / (a - b)
    else:
        return (div_diff(f_, xs_[1:]) - div_diff(f_, xs_[0:-1])) / (xs_[-1] - xs_[0])


def muller_method(function: Func, x0: Num, x1: Num, x2: Num, iterations: int,
                  get_full_result: bool = True) -> FLOAT_LIST:
    """Return the root calculated using Muller's method."""

    root_ = [x0, x1, x2]

    for _ in range(iterations):
        w = div_diff(function, [x2, x1]) + div_diff(function, [x2, x0]) - div_diff(function, [x2, x1])
        s_delta = sqrt(w**2 - 4 * function(x2) * div_diff(function, [x2, x1, x0]))
        denominators = [w + s_delta, w - s_delta]
        # Take the higher-magnitude denominator
        x3 = x2 - 2 * function(x2) / max(denominators, key=abs)
        # Advance
        x0, x1, x2 = x1, x2, x3

        root_.append(x3)

    return root_ if get_full_result else root_[-1]


def generalized_secant_method(function, x_0, x_1, tolerance: float = TOLERANCE,
                              get_full_result: bool = True) -> FLOAT_LIST:
    x_2 = x_1 - function(x_1) / div_diff(function, [x_0, x_1])

    root_ = [x_0, x_1, x_2]

    while True:
        f0 = function(root_[-1])

        f1 = div_diff(function, [x_2, x_1])
        f2 = div_diff(function, [x_2, x_1, x_0]) * (x_2 - x_1)

        f_ = f1 + f2
        x_3 = x_2 - (f0 / f_)

        x_0, x_1, x_2 = x_1, x_2, x_3
        root_.append(x_3)

        if abs(f0) < tolerance:
            break

    return root_ if get_full_result else root_[-1]


def sidi_method(function: Callable, x_0, x_1, tolerance=TOLERANCE, get_full_result: bool = True) -> FLOAT_LIST:
    return generalized_secant_method(function, x_0, x_1, tolerance, get_full_result)


def ridder_method(function: Callable, x_0: float, x_1: float, tolerance: float = TOLERANCE,
                  get_full_result: bool = True) -> Union[float, List]:
    """
    Find the root of a function using Ridder's method.

    Parameters
    ----------
    function : Callable
        The target function for root finding.
    x_0 : float
        The lower bound of the initial bracket.
    x_1 : float
        The upper bound of the initial bracket.
    tolerance : float, optional
        Tolerance for the convergence criterion. Default is 1e-6.
    get_full_result : bool, optional
        If True, returns the list of all intermediate roots. If False, returns only the final root.
        Default is True.

    Returns
    -------
    Union[float, List]
        List of roots found during the iterations if ``get_full_result`` is True, else returns the final root.
    """

    def sign_function(functions_to_evaluate, value):
        f1, f2 = functions_to_evaluate
        return value if f1 - f2 > 0 else -value if f2 - f1 > 0 else 0

    def is_sane():
        return True if f(x_0) * f(x_1) < 0 else False

    f, root_ = function, [x_0, x_1, (x_0 + x_1) / 2]

    while is_sane():
        x_mid = (x_0 + x_1) / 2

        num_ = f(x_mid) * (x_mid - x_0)
        den_ = sqrt(f(x_mid)**2 - f(x_0) * f(x_1))
        f_ = num_ / den_

        x_new = x_mid + sign_function([f(x_0), f(x_1)], f_)

        if abs(f(x_new)) < tolerance or f_ == 0:
            break

        if f(x_mid) * f(x_new) < 0:
            x_0, x_1 = x_mid, x_new
        else:
            x_0 = x_new

        root_.append(x_new)

    return root_ if get_full_result else root_[-1]


def steffensen_method(function: Callable, x_0: float, x_1: float, tolerance: float = TOLERANCE,
                      get_full_result: bool = True) -> Union[float, List]:
    """
    Find the root of a function using Steffensen method.

    Parameters
    ----------
    function : Callable
        The target function for root finding.
    x_0 : float
        The lower bound of the initial bracket.
    x_1 : float
        The upper bound of the initial bracket.
    tolerance : float, optional
        Tolerance for the convergence criterion. Default is 1e-6.
    get_full_result : bool, optional
        If True, returns the list of all intermediate roots. If False, returns only the final root.
        Default is True.

    Returns
    -------
    Union[float, List]
        List of roots found during the iterations if ``get_full_result`` is True, else returns the final root.
    """

    def is_sane():
        return True if f(x_0) * f(x_1) < 0 else False

    f, root_ = function, [x_0, x_1, (x_0 + x_1) / 2]

    solve = True if is_sane() else False

    while solve:
        f_mid = f(root_[-1])
        f_mid2 = f(root_[-1] + f_mid)

        num_ = f_mid**2
        den_ = f_mid2 - f_mid

        x_new = root_[-1] - (num_ / den_)

        root_.append(x_new)

        if f(x_new) == 0 or abs(f(x_new)) < tolerance:
            break

    return root_ if get_full_result else root_[-1]
