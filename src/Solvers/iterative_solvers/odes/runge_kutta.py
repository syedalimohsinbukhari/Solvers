"""Runge Kutta Methods for single ODEs

This module implements runge-kutta methods for single ODEs. The RK methods implemented are,

- rk2_solver:
    RK2 method for odes
- rk3_solver:
    RK3 method for odes
- rk4_solver:
    RK4 method for odes

Created on Dec 16 03:28:44 2023
"""

__all__ = ['rk2_solver', 'rk3_solver', 'rk4_solver']

from typing import Callable

import numpy as np
from custom_inherit import doc_inherit

from ... import DOC_STYLE, IFloat, NdArray2


def rk2_solver(ode: Callable, x_initial: IFloat, y_initial: IFloat, step_size: IFloat = 0.1, x_max: IFloat = 1.0,
               n_decimal: IFloat = 9) -> NdArray2:
    """
    Solve a first-order ordinary differential equation using the Runge-Kutta (RK2) method.

    Parameters
    ----------
    ode:
        The ordinary differential equation in the form dy/dx = ode(x, y).
    x_initial:
        The x value in the initial condition of the ODE.
    y_initial:
        The y value in the initial condition of the ODE.
    step_size:
        Step size for numerical integration, default is 0.1.
    x_max:
        Maximum x-value for integration, default is 1.0.
    n_decimal:
        Number of decimal places to round up to. Default is 8.

    Returns
    -------
        Tuple containing arrays of x and y values at each step.
    """

    num_steps = int((x_max - x_initial) / step_size) + 1

    x_values, y_values = np.zeros(num_steps), np.zeros(num_steps)
    x_values[0], y_values[0] = x_initial, y_initial

    for i in range(1, num_steps):
        x_i, y_i = x_values[i - 1], y_values[i - 1]

        k1 = step_size * ode(x_i, y_i)
        k2 = step_size * ode(x_i + step_size, y_i + k1)

        y_values[i] = np.round(y_i + 0.5 * (k1 + k2), n_decimal)
        x_values[i] = np.round(x_i + step_size, n_decimal)

    return x_values, y_values


@doc_inherit(rk2_solver, style=DOC_STYLE)
def rk3_solver(ode: Callable, x_initial: IFloat, y_initial: IFloat, step_size: IFloat = 0.1, x_max: IFloat = 1.0,
               n_decimal: IFloat = 9) -> NdArray2:
    """Solve a first-order ordinary differential equation using the Runge-Kutta (RK3) method."""

    num_steps = int((x_max - x_initial) / step_size) + 1

    x_values, y_values = np.zeros(num_steps), np.zeros(num_steps)
    x_values[0], y_values[0] = x_initial, y_initial

    for i in range(1, num_steps):
        x_i, y_i = x_values[i - 1], y_values[i - 1]

        k1 = step_size * ode(x_i, y_i)
        k2 = step_size * ode(x_i + step_size / 2, y_i + k1 / 2)
        k3 = step_size * ode(x_i + step_size, y_i - k1 + 2 * k2)

        y_values[i] = np.round(y_i + (1 / 6) * (k1 + 4 * k2 + k3), n_decimal)
        x_values[i] = np.round(x_i + step_size, n_decimal)

    return x_values, y_values


@doc_inherit(rk2_solver, style=DOC_STYLE)
def rk4_solver(ode: Callable, x_initial: IFloat, y_initial: IFloat, step_size: IFloat = 0.1, x_max: IFloat = 1.0,
               n_decimal: IFloat = 9) -> NdArray2:
    """Solve a first-order ordinary differential equation using the Runge-Kutta (RK4) method."""

    num_steps = int((x_max - x_initial) / step_size) + 1

    x_values, y_values = np.zeros(num_steps), np.zeros(num_steps)
    x_values[0], y_values[0] = x_initial, y_initial

    for i in range(1, num_steps):
        x_i, y_i = x_values[i - 1], y_values[i - 1]

        k1 = step_size * ode(x_i, y_i)
        k2 = step_size * ode(x_i + step_size / 2, y_i + k1 / 2)
        k3 = step_size * ode(x_i + step_size / 2, y_i + k2 / 2)
        k4 = step_size * ode(x_i + step_size, y_i + k3)

        y_values[i] = np.round(y_i + (1 / 6) * (k1 + 2 * k2 + 2 * k3 + k4), n_decimal)
        x_values[i] = np.round(x_i + step_size, n_decimal)

    return x_values, y_values
