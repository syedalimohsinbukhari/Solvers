"""Predictor-Corrector Methods.

This module provides functionality to get iterative solutions using predictor-corrector method. Currently, this module
provides,

- euler_trapezoidal: To solve first-order ODE via predictor-corrector method.

In addition, this module also provides, separately,

- euler_method: To solve a first-order ODE without predictor-corrector method.

Created on Dec 20 10:55:46 2023
"""

from custom_inherit import doc_inherit

from ... import DOC_STYLE, FList, Func, IFloat, LList, N_DECIMAL, OptIFloat
from ...__backend.core_helpers_ import linear_list, round_value_
from ...__backend.errors_ import XFindNotDefined


def trapezoidal_method(ode: Func, x_current: IFloat, y_current: IFloat, x_next: IFloat, y_next: IFloat,
                       step_size: IFloat, n_decimal: int = N_DECIMAL) -> IFloat:
    """
    Calculate the next value using the Trapezoidal Rule in numerical integration.

    Parameters
    ----------
    ode:
        The ordinary differential equation function.
    x_current:
        Current x value.
    y_current:
        Current y value.
    x_next:
        Next x value.
    y_next:
        Next y value.
    step_size:
        Step size.
    n_decimal:
        Number of decimals to round the result to (default is N_DECIMAL).

    Returns
    -------
        The next y value calculated using the Trapezoidal Rule.
    """

    def _trapezoidal():
        f_n = ode(x_current, y_current)
        f_p = ode(x_next, y_next)
        return y_current + (step_size / 2) * (f_n + f_p)

    return round_value_(_trapezoidal(), n_decimal)


def euler_method(ode: Func, x_initial: IFloat, y_initial: IFloat, step_size: IFloat = 0.1, x_find: OptIFloat = None,
                 n_decimal: int = 6, full_result: bool = False) -> LList or FList:
    """
    Solve a first-order ordinary differential equation using the Euler method.

    Parameters
    ----------
    ode : Callable
        The ordinary differential equation in the form dy/dx = ode(x, y).
    x_initial : float
        Initial x-value.
    y_initial : float
        Initial y-value.
    step_size : float, optional
        Step size for numerical integration. If not provided, it is calculated based on num_steps and x_max.
    x_find : float, optional
        The value of x at which the value of function is to be determined.
    n_decimal : int, optional
        Number of decimal places to round the result, default is 6.
    full_result : bool, optional
        If True, return a tuple of arrays (x_values, y_values). If False, return the final y-value only.

    Returns
    -------
        Result of the Euler method. If full_result is True, returns a tuple of arrays (x_values, y_values).
        If full_result is False, returns the final y-value only.
    """

    if x_find is None:
        raise XFindNotDefined('The value for x to find is not defined.')

    num_steps = int((x_find - x_initial) / step_size) + 1

    x_values, y_values = [x_initial] * num_steps, [y_initial] * num_steps

    for i in range(1, num_steps):
        y_values_ = y_values[i - 1] + step_size * ode(x_values[i - 1], y_values[i - 1])
        x_values_ = x_values[i - 1] + step_size

        y_values[i] = round_value_(y_values_, n_decimal)
        x_values[i] = round_value_(x_values_, n_decimal)

        x_initial += step_size

    return [x_values, y_values] if full_result else y_values[-1]


@doc_inherit(euler_method, style=DOC_STYLE)
def euler_trapezoidal(ode: Func, x_initial: IFloat, y_initial: IFloat, step_size: IFloat = 0.1,
                      x_find: OptIFloat = None, n_decimal: int = 6, full_result: bool = False):
    """Solve a first-order ordinary differential equation using the Euler method with trapezoidal correction."""

    if x_find is None:
        raise XFindNotDefined('The value for x to find is not defined.')

    num_steps = int((x_find - x_initial) / step_size) + 1

    x_values = linear_list(x_initial, x_initial + (num_steps * step_size), num_steps + 1, n_decimal=n_decimal)

    y_values = [0] * len(x_values)

    predictor_values, corrector_values = [0] * len(x_values), [0] * len(x_values)
    y_values[0], predictor_values[0], corrector_values[0] = [y_initial] * 3

    for i in range(1, num_steps + 1):
        predictor = euler_method(ode, x_values[i - 1], corrector_values[i - 1], step_size, n_decimal)
        predictor_values[i] = round_value_(predictor, n_decimal)

        corrector = trapezoidal_method(ode, x_values[i - 1], y_values[i - 1], x_values[i], predictor_values[i],
                                       step_size, n_decimal)
        corrector_values[i] = corrector
        y_values[i] = corrector

    return (x_values, (predictor_values, corrector_values)) if full_result else (x_values, corrector_values)
