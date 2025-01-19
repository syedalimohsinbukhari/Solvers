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

from typing import List

import numpy as np
from custom_inherit import doc_inherit

from ... import DOC_STYLE, FList, IFloat, LFunc, LList, NdArray
from ...__backend.core_helpers_ import num_steps_


def rk2_multi_ode(odes: LFunc, initial_conditions: FList, step_size: IFloat = 0.1, x_max: IFloat = 1.0,
                  n_decimal: int = 9) -> List[np.ndarray]:
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
    List[np.ndarray]
        List of arrays containing the values of each quantity at each step.
    """
    if isinstance(initial_conditions, float):
        raise TypeError("The initial conditions must be a list.")

    initial_conditions = np.array(initial_conditions) if isinstance(initial_conditions, list) else initial_conditions

    num_odes = len(odes) + 1
    num_steps = num_steps_(initial_conditions.tolist()[0], x_max, step_size)

    result = np.zeros((num_steps, num_odes))
    result[0] = initial_conditions

    for i in range(1, num_steps):
        quantities = result[i - 1]

        temp_ = _k_0(num_odes, odes, quantities, step_size, rk_type=2)
        temp_ = _k_n(odes, quantities, num_odes, step_size, temp_, index=0)

        result = _multi_rk_result(num_odes, quantities, step_size, i, temp_, result, rk_type=2, n_decimal=n_decimal)

    return [result[:, i] for i in range(num_odes)]


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

        temp_ = _k_0(num_odes, odes, quantities, step_size, rk_type=3)

        temp__ = [quantities[0] + step_size / 2] + list(quantities[1:] + temp_[0, :] / 2)
        temp_[1, :] = step_size * np.array([odes[j](*temp__) for j in range(num_odes - 1)])

        temp__ = [quantities[0] + step_size] + list(quantities[1:] - temp_[0, :] + 2 * temp_[1, :])
        temp_[2, :] = step_size * np.array([odes[j](*temp__) for j in range(num_odes - 1)])

        result = _multi_rk_result(num_odes, quantities, step_size, i, temp_, result, rk_type=3, n_decimal=n_decimal)

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

        result = _multi_rk_result(num_odes, quantities, step_size, i, temp_, result, rk_type=4, n_decimal=n_decimal)

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
    k0_[0] = step_size * np.array([ode(*quantities) for ode in odes])

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
    temp__ = np.array([quantities[0] + step_size] + list(quantities[1:] + temp_[index, :]))
    temp_[index + 1] = step_size * np.array([odes[j](*temp__) for j in range(num_odes - 1)])

    return temp_


def _multi_rk_result(num_odes: int, quantities: np.ndarray, step_size: IFloat, index: int, temp_: np.ndarray,
                     result: np.ndarray, rk_type: int, n_decimal: int) -> np.ndarray:
    """
    Returns the RK method results using vectorized operations in NumPy.

    Parameters
    ----------
    num_odes: int
        The number of ODEs.
    quantities: np.ndarray
        The array holding the base values to be worked on.
    step_size: float
        The step size for the RK method.
    index: int
        The iterative index for the RK method.
    temp_: np.ndarray
        The temporary array holding the processed values.
    result: np.ndarray
        The array holding the result values.
    rk_type: int
        The type of RK method used, 2 for RK2, 3 for RK3, and 4 for RK4.
    n_decimal: int
        The number of digits to round off to.

    Returns
    -------
    np.ndarray
        The result of RK method calculations.
    """
    result[index, 0] = np.round(quantities[0] + step_size, n_decimal)

    if rk_type == 2:
        coefficients = np.array([1 / 2, 1 / 2])
    elif rk_type == 3:
        coefficients = np.array([1 / 6, 4 / 6, 1 / 6])
    else:  # rk_type == 4
        coefficients = np.array([1 / 6, 2 / 6, 2 / 6, 1 / 6])

    for j in range(1, num_odes):
        temp__ = temp_[:, j - 1]
        result[index, j] = quantities[j] + np.sum(coefficients * temp__)
        result[index, j] = np.round(result[index, j], n_decimal)

    return result
