"""Created on Nov 02 00:25:14 2023"""

from src.interpolators.interpolation.INTERPOLATION_ import INTERPOLATION


class _BaseInterpolation(INTERPOLATION):

    def __init__(self, given_values, value_to_approximate, function=None, function_values=None):
        super().__init__(given_values, value_to_approximate, function, function_values)

    def _class_check(self):
        return 'Fwd' if self.__class__.__name__ == 'FwdInterpolation' else 'Bkw'

    def difference_table(self):
        table_limit = len(self.given_values) - 1

        difference_table = [self.function_values]

        for i in range(table_limit):
            temp_ = []
            for j, k in zip(difference_table[-1][:-1], difference_table[-1][1:]):
                temp_.append(round(k - j, 4))

            difference_table.append(temp_)

        return difference_table[1:]

    def interpolate(self):
        def find_p():
            approx, given = self.value_to_approximate, self.given_values
            num_ = approx - given[0] if self._class_check() == 'Fwd' else approx - given[-1]
            den_ = given[1] - given[0]

            return num_ / den_

        def factorial(number):
            return 1 if number == 0 else number * factorial(number - 1)

        idx_ = 0 if self._class_check() == 'Fwd' else -1
        difference_table = self.difference_table()
        initial_value = self.function_values[idx_]

        result = [initial_value]
        iter_condition = len(self.given_values) - 1

        _, p_value = find_p(), [find_p()]
        for i in range(1, iter_condition):
            _ *= find_p() - i if self._class_check() == 'Fwd' else find_p() + i
            p_value.append(_)

        for i in range(iter_condition):
            value = difference_table[-1][idx_] if i == iter_condition - 1 else difference_table[i][idx_]

            result.append((value * p_value[i]) / factorial(i + 1))

        return result, sum(result)


class FwdInterpolation(_BaseInterpolation):
    pass


class BkwInterpolation(_BaseInterpolation):
    pass


class DividedInterpolation(_BaseInterpolation):
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

        return sum(n_polynomial)
