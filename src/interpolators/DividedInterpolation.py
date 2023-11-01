"""Created on Oct 29 18:18:11 2023"""

import numpy as np

from newtonian.ERRORS_ import AtLeastOneParameterRequired


class DividedInterpolation:
    """
    A class for performing divided difference interpolation using Newton's method.

    Parameters
    ----------
    given_values : array-like
        The x-values for interpolation.
    value_to_approximate : float
        The x-value at which you want to approximate the function.
    function : callable, optional
        A function that returns the corresponding y-values for given x-values.
    function_values : array-like, optional
        Precomputed function values corresponding to the given x-values.

    Raises
    ------
    AtLeastOneParameterRequired
        If both `function` and `function_values` are missing.

    Attributes
    ----------
    given_values : array-like
        The given x-values for interpolation.
    value_to_approx : float
        The x-value at which you want to approximate the function.
    function_values : array-like
        The y-values corresponding to the given x-values.
    """

    def __init__(self, given_values, value_to_approximate, function=None, function_values=None):
        self.given_values = given_values
        self.value_to_approx = value_to_approximate

        if function is None and function_values is None:
            raise AtLeastOneParameterRequired("One of `function` or `function_values` parameter is required.")

        self.function_values = function_values if function_values else [function(value) for value in given_values]

    def difference_table(self, complete_table=False):
        """
        Calculate the divided difference table for interpolation.

        Returns
        -------
        np.ndarray
            The top row for the divided difference table.
        """

        n = len(self.function_values)
        difference_table = np.zeros((n, n))
        difference_table[:, 0] = self.function_values

        for j in range(1, n):
            for i in range(n - j):
                numerator = difference_table[i + 1, j - 1] - difference_table[i, j - 1]
                denominator = self.given_values[i + j] - self.given_values[i]
                difference_table[i, j] = numerator / denominator

        return [list(i) for i in difference_table] if complete_table else difference_table[0, :]

    def solve(self):
        """
        Perform divided difference interpolation to approximate the function.

        Returns
        -------
        float
            The interpolated value at `value_to_approx`.
        """
        n = len(self.function_values)

        product, all_products = 1, [1]
        for i in range(1, n):
            product *= self.value_to_approx - self.given_values[i - 1]
            all_products.append(product)

        difference_table = self.difference_table()
        n_polynomial = difference_table * np.array(all_products)

        return sum(n_polynomial)
