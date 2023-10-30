"""Created on Oct 19 03:46:05 2023"""

from ERRORS_ import AtLeastOneParameterRequired


class INTERPOLATION:

    def __init__(self, given_values, value_to_approximate, function=None, function_values=None):
        self.given_values = given_values
        self.value_to_approx = value_to_approximate

        if function is None and function_values is None:
            raise AtLeastOneParameterRequired("One of `function` or `function_values` parameter is required.")

        self.function = function if function else None
        self.function_values = function_values if function_values else None

    def _class(self):
        return 'Fwd' if self.__class__.__name__ == 'FwdInterpolation' else 'Bkw'

    def difference_table(self, table_limit=None):
        table_limit = len(self.given_values) - 1 if not table_limit else table_limit

        difference_table = [self.function_values]

        for i in range(table_limit):
            temp_ = []
            for j, k in zip(difference_table[-1][:-1], difference_table[-1][1:]):
                temp_.append(round(k - j, 4))

            difference_table.append(temp_)

        return difference_table[1:]

    def solve(self):
        def factorial(number):
            return 1 if number == 0 else number * factorial(number - 1)

        idx_ = 0 if self._class() == 'Fwd' else -1
        difference_table = self.difference_table()
        initial_value = self.function_values[idx_]

        result = [initial_value]
        iter_condition = len(self.given_values) - 1

        _, p_value = self._p(), [self._p()]
        for i in range(1, iter_condition):
            _ *= self._p() - i if self._class() == 'Fwd' else self._p() + i
            p_value.append(_)

        for i in range(iter_condition):
            value = difference_table[-1][idx_] if i == iter_condition - 1 else difference_table[i][idx_]

            result.append((value * p_value[i]) / factorial(i + 1))

        return result, sum(result)

    def _p(self):
        approx, given = self.value_to_approx, self.given_values
        num_ = approx - given[0] if self._class() == 'Fwd' else approx - given[-1]
        den_ = given[1] - given[0]

        return num_ / den_
