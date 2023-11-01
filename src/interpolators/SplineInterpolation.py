"""Created on Nov 01 23:03:42 2023"""

from src.interpolators.interpolation.INTERPOLATION_ import INTERPOLATION


class LinearSpline(INTERPOLATION):

    def __init__(self, given_values, value_to_approximate, function=None, function_values=None):
        super().__init__(given_values, value_to_approximate, function, function_values)

    def interpolate(self):
        data_points = len(self.given_values)
        for i in range(data_points - 1):
            if self.given_values[i] <= self.value_to_approx <= self.given_values[i + 1]:
                x0, x1 = self.given_values[i], self.given_values[i + 1]
                y0, y1 = self.function_values[i], self.function_values[i + 1]

                fraction1 = (self.value_to_approx - x1) / (x0 - x1)
                fraction1 *= y0

                fraction2 = (self.value_to_approx - x0) / (x1 - x0)
                fraction2 *= y1

                return fraction1 + fraction2
