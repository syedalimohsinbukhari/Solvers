"""Runge Kutta Methods for Multi-ODE Systems

This module implements runge-kutta methods for Multi-ODE systems. The RK methods implemented are,

- rk2_multi_ode:
    RK2 method for multi-odes
- rk3_multi_ode:
    RK3 method for multi-odes
- rk4_multi_ode:
    RK4 method for multi-odes

The module also generalizes a few RK operations, such as,

- _k_0: For RK2, Rk3, and RK4
    To compute the k0 factor for RK methods.
- _k_n: For RK2, and RK4
    To compute the kn factor for RK methods.
- _multi_rk_result: For RK2, Rk3, and RK4
    Method to implement various final results of RK method depending upon the method used.

Created on Jan 05 21:54:09 2024
"""

__all__ = ['rk2_multi_ode', 'rk3_multi_ode', 'rk4_multi_ode']

import numpy as np
from custom_inherit import doc_inherit

from ... import DOC_STYLE, FList, IFloat, LFunc, LList, N_DECIMAL, NdArray
from ...__backend.extra_ import num_steps_


def rk2_multi_ode(odes: LFunc, initial_conditions: FList, step_size: IFloat = 0.1, x_max: IFloat = 1.0,
                  n_decimal: IFloat = N_DECIMAL) -> LList:
    """
    Solve a system of first-order ordinary differential equations using the Runge-Kutta (RK2) method.

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
    num_steps = num_steps_(initial_conditions[0], x_max, step_size)

    result = np.zeros((num_steps, num_odes))
    result[0] = initial_conditions

    for i in range(1, num_steps):
        quantities = result[i - 1]

        temp_ = _k_0(num_odes, odes, quantities, step_size, 2)
        temp_ = _k_n(odes, quantities, num_odes, step_size, temp_, 0)

        result = _multi_rk_result(num_odes, quantities, step_size, i, temp_, result, 2, n_decimal)

    return [result[:, i].tolist() for i in range(num_odes)]


@doc_inherit(rk2_multi_ode, style=DOC_STYLE)
def rk3_multi_ode(odes: LFunc, initial_conditions: FList, step_size: IFloat = 0.1, x_max: IFloat = 1.0,
                  n_decimal: IFloat = 9) -> LList:
    """Solve a system of first-order ordinary differential equations using the Runge-Kutta (RK3) method."""

    num_odes = len(odes) + 1
    num_steps = num_steps_(initial_conditions[0], x_max, step_size)

    result = np.zeros((num_steps, num_odes))
    result[0] = initial_conditions

    for i in range(1, num_steps):
        quantities = result[i - 1]

        temp_ = _k_0(num_odes, odes, quantities, step_size, 3)

        temp__ = [quantities[0] + step_size / 2] + list(quantities[1:] + temp_[0, :] / 2)
        temp_[1, :] = step_size * np.array([odes[j](*temp__) for j in range(num_odes - 1)])

        temp__ = [quantities[0] + step_size] + list(quantities[1:] - temp_[0, :] + 2 * temp_[1, :])
        temp_[2, :] = step_size * np.array([odes[j](*temp__) for j in range(num_odes - 1)])

        result = _multi_rk_result(num_odes, quantities, step_size, i, temp_, result, 3, n_decimal)

    return [result[:, i].tolist() for i in range(num_odes)]


@doc_inherit(rk2_multi_ode, style=DOC_STYLE)
def rk4_multi_ode(odes: LFunc, initial_conditions: FList, step_size: IFloat = 0.1, x_max: IFloat = 1.0,
                  n_decimal: IFloat = 9) -> LList:
    """Solve a system of first-order ordinary differential equations using the Runge-Kutta (RK4) method."""

    num_odes = len(odes) + 1
    num_steps = num_steps_(initial_conditions[0], x_max, step_size)

    result = np.zeros((num_steps, num_odes))
    result[0] = initial_conditions

    for i in range(1, num_steps):
        quantities = result[i - 1]

        temp_ = _k_0(num_odes, odes, quantities, step_size, 4)

        for k in range(1, 3):
            temp__ = [quantities[0] + step_size / 2] + list(quantities[1:] + temp_[k - 1, :] / 2)
            temp_[k, :] = step_size * np.array([odes[j](*temp__) for j in range(num_odes - 1)])

        temp_ = _k_n(odes, quantities, num_odes, step_size, temp_, 2)

        result = _multi_rk_result(num_odes, quantities, step_size, i, temp_, result, 4, n_decimal)

    return [result[:, i].tolist() for i in range(num_odes)]


def _k_0(num_odes: int, odes: LFunc, quantities: NdArray, step_size: IFloat, rk_type: int) -> NdArray:
    """
    Computes the k0 parameter for RK methods.

    Parameters
    ----------
    num_odes:
        Number of ODEs.
    odes:
        A list containing the ODEs.
    quantities:
        The values used for computing the k0 factor
    step_size:
        The step size for Rk method
    rk_type:
        The type of RK method used, 2 for RK2, 3 for RK3, and 4 for RK4 method.

    Returns
    -------
        k0 parameter of RK method.
    """

    k0_ = np.zeros((rk_type, num_odes - 1))

    for j in range(num_odes - 1):
        k0_[0][j] = step_size * odes[j](*quantities)

    return k0_


@doc_inherit(_k_0, style=DOC_STYLE)
def _k_n(odes: LFunc, quantities: NdArray, num_odes: int, step_size: IFloat, temp_: NdArray, index: int) -> NdArray:
    """
    Compute the kn parameter for RK methods.

    Parameters
    ----------
    temp_:
        The temporary array holding the data
    index:
        The index of the values to pick.

    Returns
    -------
        Modified temporary array.
    """
    temp__ = [quantities[0] + step_size] + list(quantities[1:] + temp_[index, :])
    temp_[index + 1, :] = step_size * np.array([odes[j](*temp__) for j in range(num_odes - 1)])

    return temp_


def _multi_rk_result(num_odes: int, quantities: NdArray, step_size: IFloat, index: int, temp_: NdArray, result: NdArray,
                     rk_type: int, n_decimal: int) -> NdArray:
    """
    Returns the RK method results.

    Parameters
    ----------
    index:
        The iterative index for the RK method.
    temp_:
        The temporary array holding the processed values.
    quantities:
        The array holding the base values to be worked on.
    step_size:
        The step size for RK method.
    num_odes:
        The number of ODEs.
    result:
        The array holding the result values.
    rk_type:
        The type of RK method used, 2 for RK2, 3 for RK3, and 4 for RK4 method.
    n_decimal:
        The number of digits to round off to.

    Returns
    -------
        The result of RK method calculations.
    """
    for j in range(num_odes):
        if j == 0:
            result[index, j] = np.round(quantities[j] + step_size, n_decimal)
        else:
            temp__ = temp_[:, j - 1]
            if rk_type == 2:
                result[index, j] = quantities[j] + (1 / 2) * (temp__[0] + temp__[1])
            elif rk_type == 3:
                result[index, j] = quantities[j] + (1 / 6) * (temp__[0] + 4 * temp__[1] + temp__[2])
            else:
                result[index, j] = quantities[j] + (1 / 6) * (temp__[0] + 2 * np.sum(temp__[1:3]) + temp__[3])
            result[index, j] = np.round(result[index, j], n_decimal)

    return result
