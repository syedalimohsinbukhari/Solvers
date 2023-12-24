"""Created on Dec 16 03:28:44 2023"""

from typing import Callable, List, Tuple

import numpy as np

NdArray2 = [np.ndarray, np.ndarray]
NdArrayN = [np.ndarray, ...]
RK2_OUTPUT = Tuple[np.ndarray, np.ndarray]


# TODO: Multi RK2, RK3 solvers and classes


def rk2_solver(ode: Callable, x_initial: float, y_initial: float, step_size: float = 0.1,
               x_max: float = 1.0) -> RK2_OUTPUT:
    """
    Solve a first-order ordinary differential equation using the Runge-Kutta (RK2) method.

    Parameters
    ----------
    ode : Callable
        The ordinary differential equation in the form dy/dx = ode(x, y).
    x_initial : float
        Initial x-value.
    y_initial : float
        Initial y-value.
    step_size : float, optional
        Step size for numerical integration, default is 0.1.
    x_max : float, optional
        Maximum x-value for integration, default is 1.0.

    Returns
    -------
    RK2_OUTPUT
        Tuple containing arrays of x and y values at each step.
    """
    num_steps = int((x_max - x_initial) / step_size) + 1
    x_values = np.zeros(num_steps)
    y_values = np.zeros(num_steps)

    x_values[0] = x_initial
    y_values[0] = y_initial

    for i in range(1, num_steps):
        x_i, y_i = x_values[i - 1], y_values[i - 1]

        k1 = step_size * ode(x_i, y_i)
        k2 = step_size * ode(x_i + step_size, y_i + k1)

        y_values[i] = y_i + 0.5 * (k1 + k2)
        x_values[i] = x_i + step_size

    return x_values, y_values


def rk3_solver(ode: Callable, x_initial: float, y_initial: float, step_size: float = 0.1,
               x_max: float = 1.0) -> NdArray2:
    """
    Solve a first-order ordinary differential equation using the Runge-Kutta (RK3) method.

    Parameters
    ----------
    ode : Callable
        The ordinary differential equation in the form dy/dx = ode(x, y).
    x_initial : float
        Initial x-value.
    y_initial : float
        Initial y-value.
    step_size : float, optional
        Step size for numerical integration, default is 0.1.
    x_max : float, optional
        Maximum x-value for integration, default is 1.0.

    Returns
    -------
    NdArray2
        Tuple containing arrays of x and y values at each step.
    """

    num_steps = int((x_max - x_initial) / step_size) + 1
    x_values = np.zeros(num_steps)
    y_values = np.zeros(num_steps)

    x_values[0] = x_initial
    y_values[0] = y_initial

    for i in range(1, num_steps):
        x_i, y_i = x_values[i - 1], y_values[i - 1]

        k1 = step_size * ode(x_i, y_i)
        k2 = step_size * ode(x_i + step_size / 2, y_i + k1 / 2)
        k3 = step_size * ode(x_i + step_size, y_i - k1 + 2 * k2)

        y_values[i] = y_i + (1 / 6) * (k1 + 4 * k2 + k3)
        x_values[i] = x_i + step_size

    return x_values, y_values


def rk4_solver(ode: Callable, x_initial: float, y_initial: float, step_size: float = 0.1,
               x_max: float = 1.0) -> NdArray2:

    num_steps = int((x_max - x_initial) / step_size) + 1
    x_values = np.zeros(num_steps)
    y_values = np.zeros(num_steps)

    x_values[0] = x_initial
    y_values[0] = y_initial

    for i in range(1, num_steps):
        x_i, y_i = x_values[i - 1], y_values[i - 1]

        k1 = step_size * ode(x_i, y_i)
        k2 = step_size * ode(x_i + step_size / 2, y_i + k1 / 2)
        k3 = step_size * ode(x_i + step_size / 2, y_i + k2 / 2)
        k4 = step_size * ode(x_i + step_size, y_i + k3)

        y_values[i] = y_i + (1 / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
        x_values[i] = x_i + step_size

    return x_values, y_values


def rk4_multi_ode(odes: List[Callable], initial_conditions: List[float], step_size: float = 0.1, x_max: float = 1.0,
                  n_round: int = 9) -> NdArrayN:
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
    n_round : int, optional
        Number of decimal places to round the result, default is 9.

    Returns
    -------
    NdArrayN
        List of arrays containing the values of each quantity at each step.
    """

    num_odes = len(odes) + 1
    num_steps = int((x_max - initial_conditions[0]) / step_size) + 1

    result = np.zeros((num_steps, num_odes))
    result[0] = initial_conditions

    for i in range(1, num_steps):
        quantities = result[i - 1]
        m = np.zeros((4, num_odes - 1))

        for j in range(num_odes - 1):
            m[0][j] = step_size * odes[j](*quantities)

        for k in range(1, 3):
            temp = [quantities[0] + step_size / 2] + list(quantities[1:] + m[k - 1, :] / 2)
            m[k, :] = step_size * np.array([odes[j](*temp) for j in range(num_odes - 1)])

        temp = [quantities[0] + step_size] + list(quantities[1:] + m[2, :])
        m[3, :] = step_size * np.array([odes[j](*temp) for j in range(num_odes - 1)])

        for j in range(num_odes):
            if j == 0:
                result[i, j] = quantities[j] + step_size
            else:
                temp = m[:, j - 1]
                result[i, j] = quantities[j] + (1 / 6) * (temp[0] + 2 * np.sum(temp[1:3]) + temp[3])
                result[i, j] = np.round(result[i, j], n_round)

    return [result[:, i] for i in range(num_odes)]
