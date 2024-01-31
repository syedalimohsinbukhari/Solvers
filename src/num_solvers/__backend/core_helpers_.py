"""Utility Functions for Numeric Operations

This module provides a collection of utility functions for numerical operations, including rounding values and lists,
generating linearly spaced lists, determining the number of steps for numerical approximation methods, and filtering
similar values from a list.

Functions:
    - round_value_: Rounds off a given value.
    - round_list_: Maps the round function to a list.
    - linear_list: Creates a linearly spaced list of floats.
    - num_steps_: Generates the steps for numerical approximation methods.
    - filter_similar_values: Remove similar values from a list within a specified tolerance.

Exceptions:
    - AtLeastOneParameterRequired: Raised when neither ``step_size`` nor ``n_elements`` is provided in the
    ``linear_list`` function.

Constants:
    - TOLERANCE: Default tolerance level for filtering similar values. Default is 1e-8.
    - N_DECIMAL: Default number of decimal places to round off to. Default is 6.

Created on Jan 28 03:52:37 2024
"""

from .errors_ import AtLeastOneParameterRequired
from .. import FList, IFloat, LList, N_DECIMAL, OptIFloat, TOLERANCE


def round_value_(original_value: IFloat, n_decimal: int = N_DECIMAL) -> IFloat:
    """
    Rounds off a given value.

    Parameters
    ----------
    original_value:
        The value to round off.
    n_decimal:
        Number of decimal places to round off to. Defaults to 8.

    Returns
    -------
    IFloat:
        Rounded off value.
    """

    return round(original_value, n_decimal)


def round_list_(original_list: FList or LList, n_decimal: int = N_DECIMAL) -> FList or LList:
    """
    Maps the round function to a list.

    Parameters
    ----------
    original_list:
        List of values to round off.
    n_decimal:
        Number of decimal places to round off to. Defaults to N_DECIMAL.

    Returns
    -------
    FList or LList:
        List of rounded floats.
    """

    def _auxiliary_function(list_of_values):
        return list(map(lambda x: round_value_(x, n_decimal), list_of_values))

    if not isinstance(original_list[0], list):
        return _auxiliary_function(original_list)
    else:
        return [_auxiliary_function(i) for i in original_list]


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
        Number of elements in the list.
    step_size:
        Step size for the list to generate. Default is None.
    n_decimal:
        Number of digits to round off the elements of the output list.

    Returns
    -------
    FList:
        Linearly spaced list of floats between ``start`` and ``stop``.
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


# TODO: See if num_steps_ can be replaced by `linear_list`

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


def filter_similar_values(original_list: FList, tolerance: int = TOLERANCE):
    """
    Remove similar values from a list within a specified tolerance.

    Parameters
    ----------
    original_list:
        The input list containing values to be filtered.
    tolerance:
        Tolerance level for considering values as similar. Default is 1e-10.

    Returns
    -------
    FList
        A new list containing unique values from the original list based on the specified tolerance.
    """

    unique_list: list = [original_list[0]]

    for value in original_list[1:]:
        if all(abs(value - unique_value) > tolerance for unique_value in unique_list):
            unique_list.append(value)

    return unique_list
