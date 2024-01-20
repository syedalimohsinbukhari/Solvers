"""Extra functionalities

This module holds general functionalities that can be used in the entire package,

- NumSolversMatrix:
    Introduced a list-type array object for num_solvers package
- round_value_:
    Rounds the given value.
- round_list_:
    Rounds the given list.
- linear_list:
    Provides a list of linearly spaced elements from start to stop.

Created on Jan 10 00:01:13 2024
"""

__all__ = ['linear_list', 'round_value_', 'round_list_', 'round_matrix_', 'num_steps_', 'filter_similar_values']

from .errors_ import AtLeastOneParameterRequired
from .. import FList, IFloat, LList, N_DECIMAL, OptIFloat, TOLERANCE
from ..matrix_decomposition.matrix import Matrix


def round_value_(value: IFloat, n_decimal: int = 8) -> IFloat:
    """
    Rounds off a given value.

    Parameters
    ----------
    value:
        The value to round off.
    n_decimal:
        Number of decimal places to round off to.

    Returns
    -------
        Rounded off value.
    """

    return round(value, n_decimal)


def round_list_(values: FList or LList, n_decimal: int = 8) -> FList or LList:
    """
    Maps the round function to a list.

    Parameters
    ----------
    values:
        List of values to round off.
    n_decimal:
        Number of decimal places to round off to.

    Returns
    -------
        List of rounded floats.
    """

    def _auxiliary_function(list_of_values):
        return list(map(lambda x: round_value_(x, n_decimal), list_of_values))

    if not isinstance(values[0], list):
        return _auxiliary_function(values)
    else:
        return [_auxiliary_function(i) for i in values]


def round_matrix_(matrix: Matrix, n_decimal: int = N_DECIMAL) -> Matrix:
    """
    Maps the round function to a matrix.

    Parameters
    ----------
    matrix:
        Given matrix to round off.
    n_decimal:
        Number of decimal places to round off to.

    Returns
    -------
        Rounded off matrix.
    """

    matrix_ = round_list_(matrix.elements, n_decimal)
    return Matrix(matrix_)


def linear_list(start: IFloat, stop: IFloat, n_elements: OptIFloat = None, step_size: OptIFloat = None,
                n_decimal: int = 8) -> FList:
    """
    Creates a linearly spaced list of floats.

    Parameters
    ----------
    start:
        Starting value for the list.
    stop:
        Stopping value for the list.
    n_elements:
        Number of elements for the list
    step_size:
        Step size for numerical integration, default is None.
    n_decimal:
        Number of digits to round off the elements of the output list.

    Returns
    -------
        Linearly spaced list of floats between ``start`` and ``stop``
    """

    if step_size is None and n_elements is None:
        raise AtLeastOneParameterRequired('Either `step_size` or `n_elements` is required.')

    if step_size is None:
        step_size = (stop - start) / (n_elements - 1)
    elif n_elements is None:
        n_elements = (stop - start) / step_size
        n_elements += 1
        n_elements = int(n_elements)

    value = [start + (i * step_size) for i in range(n_elements)]

    return round_list_(value, n_decimal)


def num_steps_(x_initial: IFloat, x_max: IFloat, step_size: IFloat, n_decimal: int = N_DECIMAL) -> int:
    """
    Generates the steps for the numerical approximation methods.

    Parameters
    ----------
    x_initial:
        Starting value for the x variable.
    x_max:
        Maximum value for the x variable.
    step_size:
        The step-size between the starting and ending x value.
    n_decimal:
        The number of digits to round off to.

    Returns
    -------
        The number of steps required to go from ``x_initial`` to ``x_max``.
    """

    return int(round_value_(((x_max - x_initial) / step_size) + 1, n_decimal))


def filter_similar_values(original_list, tolerance=TOLERANCE):
    unique_list = [original_list[0]]

    for value in original_list[1:]:
        if all(abs(value - unique_value) > tolerance for unique_value in unique_list):
            unique_list.append(value)

    return unique_list
