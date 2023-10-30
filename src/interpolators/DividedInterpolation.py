"""Created on Oct 29 18:18:11 2023"""
import numpy as np

from newtonian.ERRORS_ import AtLeastOneParameterRequired


class DividedInterpolation:

    def __init__(self, given_values, value_to_approximate, function=None, function_values=None):
        self.given_values = given_values
        self.value_to_approx = value_to_approximate

        if function is None and function_values is None:
            raise AtLeastOneParameterRequired("One of `function` or `function_values` parameter is required.")

        self.function_values = function_values if function_values else [function(values) for values in
                                                                        given_values] if function else None

    # def difference_table(self):
    #     difference = []
    #     numerators = [self.function_values]
    #     temp_ = []
    #     temp2_ = []
    #     temp3_ = []
    #     for i in range(len(self.function_values) - 1):
    #         iter1 = numerators if i == 0 else [difference[-1]]
    #         for s, t in zip(iter1[-1], iter1[-1][1:]):
    #             temp_.append(t - s)
    #         numerators.append(temp_)
    #
    #         for j in range(len(self.given_values) - (i + 1)):
    #             temp2_.append(self.given_values[i + 1 + j] - self.given_values[j])
    #
    #         for l_, m in zip(temp_, temp2_):
    #             temp3_.append(l_ / m)
    #         difference.append(temp3_)
    #
    #         temp_ = []
    #         temp2_ = []
    #         temp3_ = []
    #
    #     return difference

    def difference_table(self):
        n = len(self.function_values)
        difference = np.zeros((n, n))
        difference[:, 0] = self.function_values

        for j in range(1, n):
            for i in range(n - j):
                num_ = difference[i + 1, j - 1] - difference[i, j - 1]
                den_ = self.given_values[i + j] - self.given_values[i]
                difference[i, j] = num_ / den_

        return list(difference[0, :])


x = [1, 1.5, 2, 3.7, 4]
y = [i**2 + 3 * i for i in x]

c = DividedInterpolation(x, 3, function_values=y)
print(c.difference_table())
