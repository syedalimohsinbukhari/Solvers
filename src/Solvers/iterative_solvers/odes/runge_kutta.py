"""Created on Dec 16 03:28:44 2023"""

from typing import Callable

import numpy as np
from custom_inherit import doc_inherit

from ... import DOC_STYLE, FList, IFloat, LFunc, LList, NdArray2


# TODO: Multi RK2, RK3 solvers and classes


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


def rk4_multi_ode(odes: LFunc, initial_conditions: FList, step_size: IFloat = 0.1, x_max: IFloat = 1.0,
                  n_decimal: IFloat = 9) -> LList:
    """
    Solve a system of first-order ordinary differential equations using the Runge-Kutta (RK4) method.

    Parameters
    ----------
    odes : List[Callable]
        List of ordinary differential equations in the form dy_i/dx = odes[i](x, y_1, y_2, ..., y_n).
    initial_conditions : List[float]
        List of initial values for each quantity (x_0, y_1, y_2, ..., y_n).
    step_size : float, optional
        Step size for numerical integration, default is 0.1.
    x_max : float, optional
        Maximum x-value for integration, default is 1.0.
    n_decimal : int, optional
        Number of decimal places to round the result, default is 9.

    Returns
    -------
        List of lists containing the values of each quantity at each step.
    """

    num_odes = len(odes) + 1
    num_steps = int((x_max - initial_conditions[0]) / step_size) + 1

    result = np.zeros((num_steps, num_odes))
    result[0] = initial_conditions

    for i in range(1, num_steps):
        quantities = result[i - 1]
        temp_ = np.zeros((4, num_odes - 1))

        for j in range(num_odes - 1):
            temp_[0][j] = step_size * odes[j](*quantities)

        for k in range(1, 3):
            temp__ = [quantities[0] + step_size / 2] + list(quantities[1:] + temp_[k - 1, :] / 2)
            temp_[k, :] = step_size * np.array([odes[j](*temp__) for j in range(num_odes - 1)])

        temp__ = [quantities[0] + step_size] + list(quantities[1:] + temp_[2, :])
        temp_[3, :] = step_size * np.array([odes[j](*temp__) for j in range(num_odes - 1)])

        for j in range(num_odes):
            if j == 0:
                result[i, j] = np.round(quantities[j] + step_size, n_decimal)
            else:
                temp__ = temp_[:, j - 1]
                result[i, j] = quantities[j] + (1 / 6) * (temp__[0] + 2 * np.sum(temp__[1:3]) + temp__[3])
                result[i, j] = np.round(result[i, j], n_decimal)

    return [result[:, i].tolist() for i in range(num_odes)]
