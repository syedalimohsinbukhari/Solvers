"""Root finding algorithms

This module provides several functions for finding the roots of given uni-variable function. An example of such
function will be

f(x) = x**2 - cos(x)

Or

f(y) = y**2 + 2*y - 1

For such uni-variable functions, the root finding methods include,

- bisection_method
- false_position_method
- regula_falsi_method
- generalized_secant_method
- muller_method
- newton_raphson_method
- ridder_method
- secant_method
- sidi_method
- steffensen_method

Created on Dec 24 12:56:39 2023
"""

__all__ = ['bisection_method', 'false_position_method', 'regula_falsi_method', 'secant_method',
           'generalized_secant_method', 'sidi_method', 'ridder_method', 'steffensen_method', 'muller_method',
           'newton_raphson_method']

from cmath import sqrt

from custom_inherit import doc_inherit

from .. import DOC_STYLE, Func, IFloat, IFloatOrFList, TOLERANCE
from ..__backend.root_finding_algorithms_ import div_diff


# TODO: modify newton_raphson method out of recursive function


def bisection_method(function: Func, x_start: IFloat, x_end: IFloat, tolerance: IFloat = TOLERANCE,
                     get_all_values: bool = False) -> IFloatOrFList:
    """
    Use the bisection method to find the root of a given function within a specified interval.

    Parameters
    ----------
    function :
        The function for which the root is being sought.
    x_start :
        The starting point of the interval.
    x_end :
        The ending point of the interval.
    tolerance :
        Tolerance for convergence. Default is 1e-10.
    get_all_values :
        If True, returns a list of all intermediate values. If False, returns only the final root.

    Returns
    -------
        The approximate root of the function or a list of all intermediate values.

    Notes
    -----
        The function assumes that the initial interval [x_start, x_end] contains a root.
    """

    f, root = function, []

    def is_sane():
        return f(x_start) * f(x_end) < 0

    while is_sane():
        x_new = (x_start + x_end) / 2

        root.append(x_new)
        f_x_new = f(x_new)

        if abs(f_x_new) < tolerance:
            break

        if f_x_new < 0:
            x_start = x_new
        else:
            x_end = x_new

    return root if get_all_values else root[-1]


@doc_inherit(bisection_method, style=DOC_STYLE)
def false_position_method(function: Func, x_start: IFloat, x_end: IFloat, tolerance: IFloat = TOLERANCE,
                          get_all_values: bool = False) -> IFloatOrFList:
    """Use the false position method to find the root of a given function within a specified interval."""
    f, root = function, []

    def is_sane():
        return f(x_start) * f(x_end) < 0

    while is_sane():
        f_x0, f_x1 = f(x_start), f(x_end)

        x_new = (x_start * f_x1 - x_end * f_x0) / (f_x1 - f_x0)
        root.append(x_new)

        f_x_new = f(x_new)

        if abs(f_x_new) < tolerance:
            break

        if f_x_new == 0:
            return root if get_all_values else root[-1]

        if f_x0 * f_x_new < 0:
            x_end = x_new
        else:
            x_start = x_new

    return root if get_all_values else root[-1]


@doc_inherit(bisection_method, style=DOC_STYLE)
def regula_falsi_method(function: Func, x_start: IFloat, x_end: IFloat, tolerance: IFloat = TOLERANCE,
                        get_all_values: bool = False) -> IFloatOrFList:
    """Use the regula-falsi method to find the root of a given function within a specified interval."""
    return false_position_method(function, x_start, x_end, tolerance, get_all_values)


@doc_inherit(bisection_method, style=DOC_STYLE)
def secant_method(function: Func, x_start: IFloat, x_end: IFloat, tolerance: IFloat = TOLERANCE,
                  get_all_values: bool = False) -> IFloatOrFList:
    """Apply the secant method to find the root of a given function within a specified interval."""

    f, root = function, []

    while True:
        x_new = x_end - f(x_end) * (x_end - x_start) / (f(x_end) - f(x_start))

        root.append(x_new)
        x_start, x_end = x_end, x_new

        f_x_new = f(x_new)
        if abs(f_x_new) < tolerance:
            break

        if f_x_new == 0:
            return root if get_all_values else root[-1]

    return root if get_all_values else root[-1]


@doc_inherit(bisection_method, style=DOC_STYLE)
def generalized_secant_method(function: Func, x_start: IFloat, x_end: IFloat, tolerance: IFloat = TOLERANCE,
                              get_full_result: bool = True) -> IFloatOrFList:
    """Use the generalized secant method to find the root of a given function within a specified interval."""

    x_2 = x_end - function(x_end) / div_diff(function, [x_start, x_end])

    root_ = [x_start, x_end, x_2]

    while True:
        f0 = function(root_[-1])

        f1 = div_diff(function, [x_2, x_end])
        f2 = div_diff(function, [x_2, x_end, x_start]) * (x_2 - x_end)

        f_ = f1 + f2
        x_3 = x_2 - (f0 / f_)

        x_start, x_end, x_2 = x_end, x_2, x_3
        root_.append(x_3)

        if abs(f0) < tolerance:
            break

    return root_ if get_full_result else root_[-1]


@doc_inherit(bisection_method, style=DOC_STYLE)
def sidi_method(function: Func, x_start: IFloat, x_end: IFloat, tolerance: IFloat = TOLERANCE,
                get_full_result: bool = True) -> IFloatOrFList:
    """Use the Sidi method to find the root of a given function within a specified interval."""
    return generalized_secant_method(function, x_start, x_end, tolerance, get_full_result)


@doc_inherit(bisection_method, style=DOC_STYLE)
def ridder_method(function: Func, x_start: IFloat, x_end: IFloat, tolerance: IFloat = TOLERANCE,
                  get_full_result: bool = True) -> IFloatOrFList:
    """Use the Ridder method to find the root of a given function within a specified interval."""

    def sign_function(functions_to_evaluate, value):
        f1, f2 = functions_to_evaluate
        return value if f1 - f2 > 0 else -value

    def is_sane():
        return f(x_start) * f(x_end) < 0

    f, root_ = function, [x_start, x_end, (x_start + x_end) / 2]

    while is_sane():
        x_mid = (x_start + x_end) / 2

        num_ = f(x_mid) * (x_mid - x_start)
        den_ = sqrt(f(x_mid)**2 - f(x_start) * f(x_end))
        f_ = num_ / den_

        x_new = x_mid + sign_function([f(x_start), f(x_end)], f_)

        if abs(f(x_new)) < tolerance or f_ == 0:
            break

        if f(x_mid) * f(x_new) < 0:
            x_start, x_end = x_mid, x_new
        else:
            x_start = x_new

        root_.append(x_new)

    return root_ if get_full_result else root_[-1]


@doc_inherit(bisection_method, style=DOC_STYLE)
def steffensen_method(function: Func, x_start: IFloat, x_end: IFloat, tolerance: IFloat = TOLERANCE,
                      get_full_result: bool = False) -> IFloatOrFList:
    """Use the Steffensen method to find the root of a given function within a specified interval."""

    def is_sane():
        return f(x_start) * f(x_end) < 0

    f, root_ = function, [x_start, x_end, (x_start + x_end) / 2]

    while is_sane():
        f_mid = f(root_[-1])
        f_mid2 = f(root_[-1] + f_mid)

        num_ = f_mid**2
        den_ = f_mid2 - f_mid

        x_new = root_[-1] - (num_ / den_)

        root_.append(x_new)

        if f(x_new) == 0 or abs(f(x_new)) < tolerance:
            break

    return root_ if get_full_result else root_[-1]


def newton_raphson_method(function: Func, derivative_of_function: Func, initial_guess: IFloat,
                          tolerance: IFloat = TOLERANCE) -> IFloat:
    """
    Find the root of a function using the Newton-Raphson method.

    Parameters
    ----------
    function : callable
        The function for which the root is being sought.
    derivative_of_function : callable
        The derivative of the target function.
    initial_guess : float
        The initial guess for the root.
    tolerance : float, optional
        Tolerance for convergence. Default is 1e-14.

    Returns
    -------
    float
        The approximate root of the function.
    """
    f, df, x_0 = function, derivative_of_function, initial_guess

    if abs(f(x_0)) < tolerance:
        return x_0
    else:
        return newton_raphson_method(f, df, x_0 - f(x_0) / df(x_0), tolerance)


def muller_method(function: Func, x_0: IFloat, x_1: IFloat, x_2: IFloat, iterations: int,
                  get_full_result: bool = True) -> IFloatOrFList:
    """Apply the muller method to find the root of a given function within a specified interval.

    Parameters
    ----------
    function:
        The function for which the root is being sought.
    x_0:
        The first initial value.
    x_1:
        The second initial value.
    x_2:
        The third initial value.
    iterations:
        Number of iterations to perform.
    get_full_result:
        If True, returns a list of all intermediate values. If False, returns only the final root.

    Returns
    -------
        The approximate root of the function or a list of all intermediate values.

    Notes
    -----
        The function assumes that the initial interval [x_start, x_end] contains a root.
    """

    root_ = [x_0, x_1, x_2]

    for _ in range(iterations):
        w = div_diff(function, [x_2, x_1]) + div_diff(function, [x_2, x_0]) - div_diff(function, [x_2, x_1])
        s_delta = sqrt(w**2 - 4 * function(x_2) * div_diff(function, [x_2, x_1, x_0]))
        denominators = [w + s_delta, w - s_delta]
        # Take the higher-magnitude denominator
        x3 = x_2 - 2 * function(x_2) / max(denominators, key=abs)
        # Advance
        x_0, x_1, x_2 = x_1, x_2, x3

        root_.append(x3)

    return root_ if get_full_result else root_[-1]
