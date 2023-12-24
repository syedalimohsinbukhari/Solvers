"""Created on Nov 02 00:25:14 2023"""

from typing import Callable, List, Optional

from ._backend.INTERPOLATION_ import BkwInterpolation, DividedInterpolation, FwdInterpolation


class InvalidMethodError(Exception):
    pass


def newton_interpolation(method: str, given_values: List[float], value_to_approximate: float,
                         function: Optional[Callable] = None, function_values: Optional[List[float]] = None,
                         use_full_table: bool = False, get_difference_table: bool = False) -> dict:
    """
    Perform Newton Interpolation using the specified method.

    Parameters
    ----------
    method: str
        Interpolation method ('forward', 'backward', or 'divided')
    given_values: List[float]
        List of x values.
    value_to_approximate: float
        The x value for which to perform interpolation.
    function: Optional[Callable]
        The function to use for interpolation.
    function_values: Optional[List[float]]
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
        And with ``get_difference_table`` flag,
            - difference_table

    Notes
    ------
        Either ``function`` or ``function_values`` must be specified.
    """

    interpolation_classes = {'forward': FwdInterpolation,
                             'backward': BkwInterpolation,
                             'divided': DividedInterpolation}

    try:
        interpolation_class = interpolation_classes[method]
    except KeyError:
        raise InvalidMethodError("Invalid interpolation method")

    interpolation = interpolation_class(given_values, value_to_approximate, function, function_values, use_full_table)

    result = interpolation.interpolate()
    if get_difference_table:
        result['difference_table'] = interpolation.difference_table()

    return result
