"""Created on Dec 16 03:28:44 2023"""

from typing import Callable, List

import numpy as np

NdArray2 = [np.ndarray, np.ndarray]
NdArrayN = [np.ndarray, ...]


# TODO: Multi RK2, RK3 solvers and classes

def rk2_solver(ode: Callable, x_0: float, y_0: float, h: float = 0.1, x_max: float = 1.0) -> NdArray2:

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

    num_odes = len(odes) + 1
    num_steps = int((x_max - initial_conditions[0]) / h) + 1

    result = np.zeros((num_steps, num_odes))
    result[0] = initial_conditions

    for i in range(1, num_steps):
        quantities = result[i - 1]
        m = np.zeros((4, num_odes - 1))

        for j in range(num_odes - 1):
            m[0][j] = h * odes[j](*quantities)

        for k in range(1, 3):
            temp = [quantities[0] + h / 2] + list(quantities[1:] + m[k - 1, :] / 2)
            m[k, :] = h * np.array([odes[j](*temp) for j in range(num_odes - 1)])

        temp = [quantities[0] + h] + list(quantities[1:] + m[2, :])
        m[3, :] = h * np.array([odes[j](*temp) for j in range(num_odes - 1)])

        for j in range(num_odes):
            if j == 0:
                result[i, j] = quantities[j] + h
            else:
                temp = m[:, j - 1]
                result[i, j] = quantities[j] + (1 / 6) * (temp[0] + 2 * np.sum(temp[1:3]) + temp[3])
                result[i, j] = np.round(result[i, j], n_round)

    return [result[:, i] for i in range(num_odes)]
