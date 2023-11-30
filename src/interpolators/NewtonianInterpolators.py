"""Created on Nov 02 00:25:14 2023"""

from src.interpolators.interpolation.INTERPOLATION_ import INTERPOLATION


class _BaseInterpolation(INTERPOLATION):

    def __init__(self, given_values, value_to_approximate, function=None, function_values=None, use_full_table=True):
        super().__init__(given_values, value_to_approximate, function, function_values, use_full_table)

    def _class_check(self):
        return 'Fwd' if self.__class__.__name__ == 'FwdInterpolation' else 'Bkw'

    def _get_index(self):
        temp = self.given_values
        temp.insert(0, self.value_to_approximate)
        temp.sort()

        temp_idx = temp.index(self.value_to_approximate)
        del temp[temp_idx]
        return temp_idx

    def difference_table(self):

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

    def interpolate(self):
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

        return result, sum(result)


class FwdInterpolation(_BaseInterpolation):
    pass


class BkwInterpolation(_BaseInterpolation):
    pass


# x = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
# val = 1.5
#
#
# def f_of_x(_x):
#     if isinstance(_x, float):
#         return _x**3 - 5 * np.sqrt(_x)
#     else:
#         return [(i**3 - 5 * np.sqrt(i)) for i in _x]
#
#
# c1 = FwdInterpolation(x, val, function_values=f_of_x(x), use_full_table=True)
# c2 = FwdInterpolation(x, val, function_values=f_of_x(x), use_full_table=False)
# ic(f_of_x(val))
# ic(c1.interpolate()[1])
# ic(c2.interpolate()[1])
#
# ic((f_of_x(val) - c2.interpolate()[1]) / f_of_x(val))


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

# x = [0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4]
# y = [0, 4, 7.94, 11.68, 14.97, 17.39, 18.25, 16.08, 0]
# val = 0.8
#
# c = DividedInterpolation(x, val, function_values=y)
# p = c.difference_table(True)
# for i in p:
#     print(i)
# print(c.interpolate())
