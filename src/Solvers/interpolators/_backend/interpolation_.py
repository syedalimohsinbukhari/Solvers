"""Created on Oct 19 03:46:05 2023"""

from custom_inherit import doc_inherit

from . import ERRORS_
from ... import DOC_STYLE, FList, IFloat, LList, OptFunc, OptList


# TODO: Add module level docstrings


class Interpolation:
    """
    Base interpolation class.

    Methods
    -------
    difference_table:
        To calculate the difference table for the interpolation method.
    interpolate:
        Interpolates the ``value_to_approximate``.
    """

    def __init__(self, given_values: FList, value_to_approximate: IFloat, function: OptFunc = None,
                 function_values: OptList = None, n_decimal: int = 8):
        """
        Initialize the INTERPOLATION object.

        Parameters
        ----------
        given_values:
            List of x-values.
        value_to_approximate:
            The x-value to interpolate the corresponding y-value.
        function:
            The function to use for calculating y-values if ``function_values`` is not provided.
        function_values:
            List of y-values corresponding to ``given_values``. If not provided, ``function`` will be used to calculate
            these values.
        n_decimal:
            Number of digits to round off to.

        Raises
        ------
        ERRORS_.AtLeastOneParameterRequired
            If neither ``function`` nor ``function_values`` is provided.
        """

        self.given_values = given_values
        self.value_to_approximate = value_to_approximate

        if function is None and function_values is None:
            raise ERRORS_.AtLeastOneParameterRequired("One of `function` or `function_values` parameter is required.")

        self.function_values = function_values if function_values else [function(value) for value in given_values]
        self.step_size = given_values[1] - given_values[0]
        self.round_ = n_decimal

    def _class_check(self) -> str:
        """
        Check the class type.

        Returns
        -------
            'Fwd' if the class is FwdInterpolation, 'Bkw' if the class is BkwInterpolation.
        """

        return 'Fwd' if self.__class__.__name__ == 'FwdInterpolation' else 'Bkw'

    def _get_index(self) -> int:
        """
        Get the index of the value to approximate in the sorted given values.

        Returns
        -------
            The index of the value to approximate in the sorted given values.
        """

        temp = self.given_values
        temp.insert(0, self.value_to_approximate)
        temp.sort()

        temp_idx = temp.index(self.value_to_approximate)
        del temp[temp_idx]

        return temp_idx

    def difference_table(self) -> FList:
        """
        Calculate the divided difference table.

        Returns
        -------
            The divided difference table.
        """

        idx_ = self._get_index()

        table_limit = len(self.given_values) - 1

        difference_table = [self.function_values]

        for i in range(table_limit):
            temp_ = []
            for j, k in zip(difference_table[-1][:-1], difference_table[-1][1:]):
                temp_.append(round(k - j, self.round_))

            difference_table.append(temp_)

        if self._class_check() == 'Fwd':
            if idx_ - 1 != 0:
                reduced_index = idx_ - 1
                reduced_table = difference_table[:-reduced_index]
                d_table = [i[reduced_index] for i in reduced_table]
            else:
                d_table = difference_table
                d_table = [i[0] for i in d_table]
        else:
            d_table = [v[idx_ - i - 1] for i, v in enumerate(difference_table[:idx_])]

        return d_table

    def interpolate(self) -> IFloat:
        """
        Interpolate the y-value corresponding to the given x-value.

        Returns
        -------
            Dictionary with ``step_values`` (list of interpolated values for each step) and ``result``
            (sum of the interpolated values).
        """

        def find_p():
            return (approx - given[index - 1]) / step_size

        def factorial(number):
            return 1 if number == 0 else number * factorial(number - 1)

        approx, index, step_size = self.value_to_approximate, self._get_index(), self.step_size
        given, func_values, round_to, class_ = self.given_values, self.function_values, self.round_, self._class_check()

        difference_table = self.difference_table()[1:]
        initial_value = func_values[index - 1]

        result, iter_condition = [initial_value], len(difference_table)

        _, p_value = find_p(), [find_p()]

        for i in range(1, iter_condition):
            _ *= (find_p() - i) if class_ == 'Fwd' else find_p() + i
            p_value.append(_)

        result = [(difference_table[i] * p_value[i]) / factorial(i + 1) for i in range(iter_condition)]
        result.insert(0, initial_value)
        result = sum(list(map(lambda x: round(x, round_to), result)))

        return result


@doc_inherit(Interpolation, style=DOC_STYLE)
class FwdInterpolation(Interpolation):
    """Forward Interpolation class specializes in implementing the forward difference method."""


@doc_inherit(Interpolation, style=DOC_STYLE)
class BkwInterpolation(Interpolation):
    """Backward Interpolation class specializes in implementing the backward difference method."""


@doc_inherit(Interpolation.__doc__, style=DOC_STYLE)
class DividedInterpolation(Interpolation):
    """Divided Interpolation class specializes in implementing the divided difference method."""

    @doc_inherit(Interpolation.__init__, style=DOC_STYLE)
    def __init__(self, given_values: FList, value_to_approximate: IFloat, function: OptFunc = None,
                 function_values: OptList = None, n_decimal: int = 8):
        super().__init__(given_values, value_to_approximate, function, function_values, n_decimal)

    def difference_table(self) -> FList:
        given, func_values = self.given_values, self.function_values
        n_values = len(func_values)

        difference_table: LList = [[0] * n_values for _ in range(n_values)]

        for i in range(n_values):
            difference_table[i][0] = func_values[i]

        for j in range(1, n_values):
            for i in range(n_values - j):
                numerator = difference_table[i + 1][j - 1] - difference_table[i][j - 1]
                denominator = given[i + j] - given[i]
                difference_table[i][j] = numerator / denominator

        return difference_table[0]

    def interpolate(self) -> IFloat:
        approx, given, func_values = self.value_to_approximate, self.given_values, self.function_values
        n_values, round_to = len(func_values), self.round_

        product, all_products = 1, [1]
        for i in range(1, n_values):
            product *= approx - given[i - 1]
            all_products.append(product)

        difference_table = self.difference_table()
        result = round(sum(i * j for i, j in zip(difference_table, all_products)), round_to)

        return result
