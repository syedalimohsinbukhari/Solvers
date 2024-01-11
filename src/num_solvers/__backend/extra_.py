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

__all__ = ['linear_list', 'round_value_', 'round_list_', 'NumSolversMatrix']

from .errors_ import AtLeastOneParameterRequired, InconsistentDimensions, IndexCanNotBeSlice
from .. import FList, IFloat, LList, OptIFloat


class NumSolversMatrix:
    """Class for NumSolversMatrix object."""

    def __init__(self, values: LList or FList):
        """
        Initializer class for NumSolverArray.

        Parameters
        ----------
        values:
            Values to transform into a NumSolverArray object.
        """

        self.values = values
        self.__sanity_check()

    def __sanity_check(self):
        """Performs a sanity check on NumSolversMatrix object."""

        if isinstance(self.values[0], list):
            len1 = len(self.values[0])
            if all(len(x) == len1 for x in self.values[1:]):
                pass
            else:
                raise InconsistentDimensions('The members of NumSolversMatrix are not of same length.')

    def __repr__(self):
        vals = self.values
        len_max = len(str(max(max(vals)))) if isinstance(vals[0], list) else len(str(max(vals)))
        len_cond = len_max + 5 if len_max % 2 == 0 else len_max + 4

        out = '['
        for i, v in enumerate(vals):
            if isinstance(v, (int, float)):
                out += f'{v:^{len_cond}.{len_cond - 4}f}' if i + 1 == len(vals) else f'{v:^{len_cond}.{len_cond - 4}f}|'
            else:
                for k in v:
                    out += f'{k:^{len_cond}.5f}' if i + 1 == len(v) else f'{k:^{len_cond}.5f}|'
                if v != vals[-1]:
                    out += ']\n['
        out += ']'

        return out

    def __getitem__(self, item):
        return self.values[item][0] if isinstance(item, slice) else self.values[item]

    def __setitem__(self, index, value):
        if isinstance(self.values[index], list):
            raise IndexCanNotBeSlice()
        self.values[index] = value


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
