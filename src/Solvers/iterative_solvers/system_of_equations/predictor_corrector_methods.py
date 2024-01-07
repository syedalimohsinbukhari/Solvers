"""Created on Dec 20 10:55:46 2023"""

from typing import Callable, Optional, Tuple, Union

import numpy as np

from ...__backend.predictor_corrector_ import get_x, trapezoidal_rule

FLOAT = Union[float, np.float32]
EULER_OUTPUT = Union[Tuple[np.ndarray, np.ndarray], float]
EULER_TRAPEZOIDAL_OUTPUT = Union[Tuple[np.ndarray, np.ndarray], Tuple[np.ndarray, Tuple[np.ndarray, np.ndarray]]]


def euler_method(ode: Callable, x_initial: FLOAT, y_initial: FLOAT, step_size: Optional[float] = 0.1,
                 num_steps: Optional[int] = 1, x_max: Optional[float] = None, n_decimal: int = 6,
                 full_result: bool = False) -> EULER_OUTPUT:
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
    num_steps : int, optional
        Number of steps for numerical integration. Default is 1.
    x_max : float, optional
        Maximum x-value for integration. If not provided, it is calculated based on x_initial, step_size, and num_steps.
    n_decimal : int, optional
        Number of decimal places to round the result, default is 6.
    full_result : bool, optional
        If True, return a tuple of arrays (x_values, y_values). If False, return the final y-value only.

    Returns
    -------
    EULER_OUTPUT
        Result of the Euler method. If full_result is True, returns a tuple of arrays (x_values, y_values).
        If full_result is False, returns the final y-value only.
    """
    step_size = (x_max - x_initial) / num_steps if step_size is None else step_size
    x_max = x_initial + num_steps * step_size if x_max is None else x_max

    y_values, x_values = [y_initial], [x_initial]

    while x_initial < x_max:
        y_values.append(y_values[-1] + step_size * ode(x_values[-1], y_values[-1]))
        x_values.append(x_values[-1] + step_size)

        x_initial = np.round(x_initial + step_size, n_decimal)

    return (np.array(x_values), np.array(y_values)) if full_result else y_values[-1]


def euler_trapezoidal(ode: Callable, x_initial: float, y_initial: float, step_size: float, num_steps: int,
                      n_decimal: int = 6, get_predictor: bool = False) -> EULER_TRAPEZOIDAL_OUTPUT:
    """
    Solve a first-order ordinary differential equation using the Euler method with trapezoidal correction.

    Parameters
    ----------
    ode : Callable
        The ordinary differential equation in the form dy/dx = ode(x, y).
    x_initial : float
        Initial x-value.
    y_initial : float
        Initial y-value.
    step_size : float
        Step size for numerical integration.
    num_steps : int
        Number of steps for numerical integration.
    n_decimal : int, optional
        Number of decimal places to round the result, default is 6.
    get_predictor : bool, optional
        If True, return both the corrector values and the predictor-corrector pair.
        If False, return only the corrector values.

    Returns
    -------
    EULER_TRAPEZOIDAL_OUTPUT
        Result of the Euler method with trapezoidal correction.

        - If ``get_predictor`` is True, returns (x_values, (predictor_values, corrector_values)).
        - If ``get_predictor`` is False, returns (x_values, corrector_values).
    """
    x_values = get_x(x_initial, num_steps, step_size, n_decimal)
    y_values = np.zeros_like(x_values)

    predictor_values, corrector_values = np.zeros_like(x_values), np.zeros_like(x_values)
    y_values[0], predictor_values[0], corrector_values[0] = [y_initial] * 3

    for i in range(1, num_steps + 1):
        predictor = euler_method(ode, x_values[i - 1], corrector_values[i - 1], step_size, n_decimal)
        predictor_values[i] = np.round(predictor, n_decimal)

        corrector = trapezoidal_rule(ode, x_values[i - 1], y_values[i - 1], x_values[i], predictor_values[i],
                                     step_size, n_decimal)
        corrector_values[i] = corrector
        y_values[i] = corrector

    return (x_values, (predictor_values, corrector_values)) if get_predictor else (x_values, corrector_values)
