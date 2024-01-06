"""Created on Nov 01 16:47:03 2023"""

from ._backend.interpolation_ import Interpolation


class LagrangeInterpolation(Interpolation):

    def __init__(self, given_values, value_to_approximate, function=None, function_values=None):
        super().__init__(given_values, value_to_approximate, function, function_values)

    def interpolate(self, give_values=False):
        answer = []
        given = self.given_values
        for i in range(len(given)):
            modified_array = given[:i] + given[i + 1:]
            numerator = 1
            denominator = 1
            for value in modified_array:
                numerator *= (self.value_to_approximate - value)
                denominator *= (given[i] - value)
            answer.append((numerator / denominator) * self.function_values[i])

        result = sum(answer)

        return (result, answer) if give_values else result
