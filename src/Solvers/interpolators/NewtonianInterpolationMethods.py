"""Created on Nov 02 00:25:14 2023"""

from typing import Callable, Optional

from . import L_FLOAT
from ._backend.INTERPOLATION_ import BkwInterpolation, DividedInterpolation, FwdInterpolation, get_result


def forward_interpolation(given_values: L_FLOAT, value_to_approximate: float, function: Optional[Callable] = None,
                          function_values: Optional[L_FLOAT] = None, use_full_table: bool = False,
                          get_difference_table: bool = False):
    """
    Perform forward Newton Interpolation.

    Parameters
    ----------
    given_values: L_FLOAT
        List of x values.
    value_to_approximate: float
        The x value for which to perform interpolation.
    function: Optional[Callable]
        The function to use for interpolation.
    function_values: Optional[L_FLOAT]
        List of function values corresponding to given_values.
    use_full_table: bool
        Flag indicating whether to use the full difference table.
    get_difference_table: bool
        Flag indicating whether to return the difference table.

    Returns
    -------
    dict
        A dictionary containing
            - step_values
            - results
        Additionally, with ``get_difference_table`` flag,
            - difference_table

    Notes
    ------
        Either ``function`` or ``function_values`` must be specified.
    """
    interpolation = FwdInterpolation(given_values, value_to_approximate, function, function_values, use_full_table)
    return get_result(interpolation, get_difference_table)


def backward_interpolation(given_values: L_FLOAT, value_to_approximate: float, function: Optional[Callable] = None,
                           function_values: Optional[L_FLOAT] = None, use_full_table: bool = False,
                           get_difference_table: bool = False):

    interpolation = BkwInterpolation(given_values, value_to_approximate, function, function_values, use_full_table)
    return get_result(interpolation, get_difference_table)


def divided_interpolation(given_values: L_FLOAT, value_to_approximate: float, function: Optional[Callable] = None,
                          function_values: Optional[L_FLOAT] = None, use_full_table: bool = False,
                          get_difference_table: bool = False):

    interpolation = DividedInterpolation(given_values, value_to_approximate, function, function_values, use_full_table)
    return get_result(interpolation, get_difference_table)
