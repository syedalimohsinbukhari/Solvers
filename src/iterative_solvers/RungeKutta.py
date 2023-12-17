"""Created on Dec 16 03:28:44 2023"""
from typing import Callable, List

import matplotlib.pyplot as plt
import numpy as np

NdArray2 = [np.ndarray, np.ndarray]
NdArrayN = [np.ndarray, ...]


def rk2_solver(ode: Callable, x_0: float, y_0: float, h: float = 0.1, x_max: float = 1.0) -> NdArray2:
    """
    Solve a differential equation using the RK2 method.

    Parameters
    ----------
    ode: Callable
        The ordinary differential equation function.
    x_0: float
        Initial x value.
    y_0: float
        Initial y value.
    h: float, optional
        Step size. Default is 0.1.
    x_max: float, optional
        Maximum x value. Default is 1.0.

    Returns
    -------
    object: NdArray
        The x and y values over the iteration as two separate arrays
    """

    num_steps = int((x_max - x_0) / h) + 1
    x_n = np.zeros(num_steps)
    y_n = np.zeros(num_steps)

    x_n[0] = x_0
    y_n[0] = y_0

    for i in range(1, num_steps):
        x_i, y_i = x_n[i - 1], y_n[i - 1]

        k1 = h * ode(x_i, y_i)
        k2 = h * ode(x_i + h, y_i + k1)

        y_n[i] = y_i + 0.5 * (k1 + k2)
        x_n[i] = x_i + h

    return x_n, y_n


def rk3_solver(ode: Callable, x_0: float, y_0: float, h: float = 0.1, x_max: float = 1.0) -> NdArray2:
    """
    Solve a differential equation using the RK3 method.

    Parameters
    ----------
    ode : Callable
        The ordinary differential equation function.
    x_0 : float
        Initial x value.
    y_0 : float
        Initial y value.
    h : float, optional
        Step size. Default is 0.1.
    x_max : float, optional
        Maximum x value. Default is 1.0.

    Returns
    -------
    NdArray2
        Arrays containing x and y values over the iteration.
    """

    num_steps = int((x_max - x_0) / h) + 1
    x_n = np.zeros(num_steps)
    y_n = np.zeros(num_steps)

    x_n[0] = x_0
    y_n[0] = y_0

    for i in range(1, num_steps):
        x_i, y_i = x_n[i - 1], y_n[i - 1]

        k1 = h * ode(x_i, y_i)
        k2 = h * ode(x_i + h / 2, y_i + k1 / 2)
        k3 = h * ode(x_i + h, y_i - k1 + 2 * k2)

        y_n[i] = y_i + (1 / 6) * (k1 + 4 * k2 + k3)
        x_n[i] = x_i + h

    return x_n, y_n


def rk4_solver(ode: Callable, x_0: float, y_0: float, h: float = 0.1, x_max: float = 1.0) -> NdArray2:
    """
    Solve a differential equation using the RK4 method.

    Parameters
    ----------
    ode : Callable
        The ordinary differential equation function.
    x_0 : float
        Initial x value.
    y_0 : float
        Initial y value.
    h : float, optional
        Step size. Default is 0.1.
    x_max : float, optional
        Maximum x value. Default is 1.0.

    Returns
    -------
    NdArray2
        Arrays containing x and y values over the iteration.
    """

    num_steps = int((x_max - x_0) / h) + 1
    x_n = np.zeros(num_steps)
    y_n = np.zeros(num_steps)

    x_n[0] = x_0
    y_n[0] = y_0

    for i in range(1, num_steps):
        x_i, y_i = x_n[i - 1], y_n[i - 1]

        k1 = h * ode(x_i, y_i)
        k2 = h * ode(x_i + h / 2, y_i + k1 / 2)
        k3 = h * ode(x_i + h / 2, y_i + k2 / 2)
        k4 = h * ode(x_i + h, y_i + k3)

        y_n[i] = y_i + (1 / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
        x_n[i] = x_i + h

    print(x_n, y_n)

    return x_n, y_n


def rk4_multi_ode(odes: List[Callable], initial_conditions: List[float], h: float = 0.1, x_max: float = 1.0,
                  n_round: int = 9) -> NdArrayN:
    """
    Solve a system of ordinary differential equations using the RK4 method.

    Parameters
    ----------
    n_round
    odes : List[Callable]
        List of ordinary differential equation functions. The variables in each ODE must be in the same order.
    initial_conditions : List[float]
        List of initial conditions for each variable in the system.
    h : float, optional
        Step size. Default is 0.1.
    x_max : float, optional
        Maximum x value. Default is 1.0.
    n_round: int, optional
        Number of decimal places to round off to. Default is 9.

    Returns
    -------
    NdArray2:
        Arrays containing x and y values over the iteration for each variable in the system.
    """

    num_odes = len(odes) + 1
    num_steps = int((x_max - initial_conditions[0]) / h) + 1

    result = np.zeros((num_steps, num_odes))
    result[0] = initial_conditions

    for i in range(1, num_steps):
        quantities = result[i - 1]
        m = np.zeros((4, num_odes - 1))

        for j in range(num_odes - 1):
            qty1 = quantities
            m[0][j] = h * odes[j](*qty1)

        for k in range(1, 3):
            qty = [quantities[0] + h / 2] + list(quantities[1:] + m[k - 1, :] / 2)
            m[k, :] = h * np.array([odes[j](*qty) for j in range(num_odes - 1)])

        qty = [quantities[0] + h] + list(quantities[1:] + m[2, :])
        m[3, :] = h * np.array([odes[j](*qty) for j in range(num_odes - 1)])

        for j in range(num_odes):
            if j == 0:
                result[i, j] = quantities[j] + h
            else:
                qty5 = m[:, j - 1]
                result[i, j] = quantities[j] + (1 / 6) * (qty5[0] + 2 * np.sum(qty5[1:3]) + qty5[3])
                result[i, j] = np.round(result[i, j], n_round)

    return [result[:, i] for i in range(num_odes)]


#######################################################################################################################
# Example
#######################################################################################################################

# def ode1(t, x1, x2, x3, x4):
#     return x2
#
#
# def ode2(t, x1, x2, x3, x4):
#     return x3
#
#
# def ode3(t, x1, x2, x3, x4):
#     return x4
#
#
# def ode4(t, x1, x2, x3, x4):
#     return -8 * x1 + np.sin(t) * x2 - 3 * x3 + t**2
#
#
# c = rk4_multi_ode([ode1, ode2, ode3, ode4], initial_conditions=[0, 1, 2, 3, 4], h=0.01, x_max=0.5)
# plt.plot(c[0], c[1], 'r-')
# plt.plot(c[0], c[2], 'g-')
# plt.plot(c[0], c[3], 'b-')
# plt.plot(c[0], c[4], 'c-')
# plt.show()
