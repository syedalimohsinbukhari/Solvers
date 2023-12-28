"""Created on Oct 19 03:46:05 2023"""

from typing import List, Optional

from . import ERRORS_


# TODO: Correct the divided interpolation method according to fwd and bkw
# TODO: Remove ``use_full_table`` option and corresponding code


class INTERPOLATION:

    def __init__(self, given_values: List[float], value_to_approximate: float,
                 function: Optional[callable] = None, function_values: Optional[List[float]] = None,
                 use_full_table: bool = True):
        """
        Initialize the INTERPOLATION object.

        Parameters
        ----------
        given_values : List[float]
            List of x-values.
        value_to_approximate : float
            The x-value to interpolate the corresponding y-value.
        function : callable, optional
            The function to use for calculating y-values if ``function_values`` is not provided.
        function_values : List[float], optional
            List of y-values corresponding to ``given_values``. If not provided, ``function`` will be used to calculate
            these values.
        use_full_table : bool, optional
            If True, use the full difference table; if False, use only the necessary parts for interpolation.

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
        self.use_full_table = use_full_table

    def difference_table(self):
        """
        Calculate the difference table.
        """
        pass

    def interpolate(self):
        """
        Interpolate the y-value corresponding to the given x-value.
        """
        pass


class BaseInterpolation(INTERPOLATION):

    def __init__(self, given_values, value_to_approximate, function=None, function_values=None, use_full_table=True):
        """
        Initialize the BaseInterpolation object.

        Parameters
        ----------
        given_values : List[float]
            List of x-values.
        value_to_approximate : float
            The x-value to interpolate the corresponding y-value.
        function : callable, optional
            The function to use for calculating y-values if `function_values` is not provided.
        function_values : List[float], optional
            List of y-values corresponding to `given_values`. If not provided, `function` will be used to calculate
            these values.
        use_full_table : bool, optional
            If True, use the full difference table; if False, use only the necessary parts for interpolation.
        """
        super().__init__(given_values, value_to_approximate, function, function_values, use_full_table)

    def _class_check(self) -> str:
        """
        Check the class type.

        Returns
        -------
        str
            'Fwd' if the class is FwdInterpolation, 'Bkw' if the class is BkwInterpolation.
        """
        return 'Fwd' if self.__class__.__name__ == 'FwdInterpolation' else 'Bkw'

    def _get_index(self) -> int:
        """
        Get the index of the value to approximate in the sorted given values.

        Returns
        -------
        int
            The index of the value to approximate in the sorted given values.
        """
        temp = self.given_values
        temp.insert(0, self.value_to_approximate)
        temp.sort()

        temp_idx = temp.index(self.value_to_approximate)
        del temp[temp_idx]
        return temp_idx

    def difference_table(self) -> List[List[float]]:
        """
        Calculate the divided difference table.

        Returns
        -------
        List[float]
            The divided difference table.
        """
        idx_ = self._get_index()

        table_limit = len(self.given_values) - 1

        difference_table = [self.function_values]

        for i in range(table_limit):
            temp_ = []
            for j, k in zip(difference_table[-1][:-1], difference_table[-1][1:]):
                temp_.append(round(k - j, 8))

            difference_table.append(temp_)

        if not self.use_full_table:
            if self._class_check() == 'Fwd':
                if idx_ - 1 != 0:
                    reduced_index = idx_ - 1
                    reduced_table = difference_table[1:-reduced_index]
                    t = [i[reduced_index] for i in reduced_table]
                else:
                    t = difference_table[1:]
                    t = [i[0] for i in t]
            else:
                t = [v[idx_ - i - 3] for i, v in enumerate(difference_table[1:])]

        else:
            t = difference_table[1:]

        return t

    def interpolate(self) -> dict:
        """
        Interpolate the y-value corresponding to the given x-value.

        Returns
        -------
        dict
            Dictionary with ``step_values`` (list of interpolated values for each step) and ``result``
            (sum of the interpolated values).
        """

        def find_p():
            if self.use_full_table:
                approx, given = self.value_to_approximate, self.given_values
                num_ = approx - given[0] if self._class_check() == 'Fwd' else approx - given[-1]
                den_ = given[1] - given[0]

                return num_ / den_
            else:
                return self.value_to_approximate - self.given_values[self._get_index() - 1]

        def factorial(number):
            return 1 if number == 0 else number * factorial(number - 1)

        idx_ = (0 if self._class_check() == 'Fwd' else -1) if self.use_full_table else (
                self._get_index() - 1 if self._class_check() == 'Fwd' else self._get_index() - 1)

        difference_table = self.difference_table()
        initial_value = self.function_values[idx_]

        result = [initial_value]
        iter_condition = len(self.given_values) - 1 if self.use_full_table else len(self.difference_table())

        _, p_value = find_p(), [find_p()]

        for i in range(1, iter_condition):
            _ *= find_p() - i if self._class_check() == 'Fwd' else find_p() + i
            p_value.append(_)

        for i in range(iter_condition):
            if self.use_full_table:
                value = difference_table[-1][idx_] if i == iter_condition - 1 else difference_table[i][idx_]

                result.append((value * p_value[i]) / factorial(i + 1))
            else:
                result.append((difference_table[i] * p_value[i]) / factorial(i + 1))

        return {'step_values': result, 'result': sum(result)}


class FwdInterpolation(BaseInterpolation):
    """
    Forward Interpolation class, inheriting from BaseInterpolation.

    This class specializes in implementing the forward interpolation method.
    It utilizes the BaseInterpolation class for shared functionality.

    Parameters and methods are inherited from BaseInterpolation.
    """
    pass


class BkwInterpolation(BaseInterpolation):
    """
    Backward Interpolation class, inheriting from BaseInterpolation.

    This class specializes in implementing the backward interpolation method.
    It utilizes the BaseInterpolation class for shared functionality.

    Parameters and methods are inherited from BaseInterpolation.
    """
    pass


class DividedInterpolation(BaseInterpolation):
    def difference_table(self, complete_table=False):

        n = len(self.function_values)
        difference_table = [[0] * n for _ in range(n)]

        for i in range(n):
            difference_table[i][0] = self.function_values[i]

        for j in range(1, n):
            for i in range(n - j):
                numerator = difference_table[i + 1][j - 1] - difference_table[i][j - 1]
                denominator = self.given_values[i + j] - self.given_values[i]
                difference_table[i][j] = numerator / denominator

        if complete_table:
            return difference_table
        else:
            return difference_table[0]

    def interpolate(self):
        n = len(self.function_values)

        product, all_products = 1, [1]
        for i in range(1, n):
            product *= self.value_to_approximate - self.given_values[i - 1]
            all_products.append(product)

        difference_table = self.difference_table()
        n_polynomial = [i * j for i, j in zip(difference_table, all_products)]

        return {'step_values': n_polynomial, 'result': sum(n_polynomial)}


def get_result(interpolation, difference_table):
    result = interpolation.interpolate()
    if difference_table:
        result['difference_table'] = interpolation.difference_table()

    return result
