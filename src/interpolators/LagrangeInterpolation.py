"""Created on Nov 01 16:47:03 2023"""

from newtonian.ERRORS_ import AtLeastOneParameterRequired


class LagrangeInterpolation:

    def __init__(self, given_values, value_to_approximate, function=None, function_values=None):
        self.given_values = given_values
        self.value_to_approx = value_to_approximate

        if function is None and function_values is None:
            raise AtLeastOneParameterRequired("One of `function` or `function_values` parameter is required.")

        self.function_values = function_values if function_values else [function(x) for x in given_values]

    # @staticmethod
    # def round_to_significant_digits(_value, threshold=1e-8):
    #     if _value == 0.0:
    #         return 0.0
    #
    #     rounded_x = threshold * round(_value / threshold)
    #     return rounded_x

    def solve(self, give_values=False):
        answer = []
        given = self.given_values
        for i in range(len(given)):
            modified_array = given[:i] + given[i + 1:]
            numerator = 1
            denominator = 1
            for value in modified_array:
                numerator *= (self.value_to_approx - value)
                denominator *= (given[i] - value)
            answer.append((numerator / denominator) * self.function_values[i])

        result = sum(answer)

        return (result, answer) if give_values else result
