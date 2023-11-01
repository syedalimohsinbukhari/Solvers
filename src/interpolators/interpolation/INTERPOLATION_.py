"""Created on Oct 19 03:46:05 2023"""

from .ERRORS_ import AtLeastOneParameterRequired


class INTERPOLATION:

    def __init__(self, given_values, value_to_approximate, function=None, function_values=None):
        self.given_values = given_values
        self.value_to_approx = value_to_approximate

        if function is None and function_values is None:
            raise AtLeastOneParameterRequired("One of `function` or `function_values` parameter is required.")

        self.function_values = function_values if function_values else [function(value) for value in given_values]

    def difference_table(self):
        pass

    def interpolate(self):
        pass
