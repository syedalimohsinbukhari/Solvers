"""Created on Dec 24 12:56:39 2023"""

from cmath import sqrt
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


def mullers_method(function: Func, x0: Num, x1: Num, x2: Num, iterations: int,
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
