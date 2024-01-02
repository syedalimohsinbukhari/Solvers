"""Created on Dec 28 17:28:17 2023"""
from typing import Callable, List, Union

import numpy as np

from .INTERPOLATION_ import INTERPOLATION


class SPLINE(INTERPOLATION):

    def __init__(self, given_values, value_to_approximate, function=None, function_values=None):
        super().__init__(given_values, value_to_approximate, function, function_values)

    def _solve(self):
        pass

    def show_splines(self):
        pass

    @staticmethod
    def x_values(start, end, n=1_000):
        return np.linspace(start, end, n)

    @staticmethod
    def get_splines(x_val, y_val=None, point=None, spline_type='linear'):
        if spline_type == 'linear':
            num1, num2 = (point - x_val[1]) * y_val[0], (point - x_val[0]) * y_val[1]
            denominator = x_val[0] - x_val[1]

            return (num1 - num2) / denominator
        elif spline_type == 'quad':
            return y_val[0] * x_val**2 + y_val[1] * x_val + y_val[2]
        elif spline_type == 'cubic':
            return y_val[0] * x_val**3 + y_val[1] * x_val**2 + y_val[2] * x_val + y_val[3]


def get_solution(given_values, value_to_approximate, solution):
    given_values.insert(0, value_to_approximate)
    given_values.sort()
    idx_ = given_values.index(value_to_approximate) - 1
    del given_values[idx_ + 1]

    return solution[idx_]


class LinearSpline(SPLINE):

    def __init__(self, given_values: List[float], value_to_approximate: float, function: Union[None, Callable] = None,
                 function_values: Union[None, List[float]] = None) -> None:
        """
        Initialize LinearSpline object.

        Parameters
        ----------
        given_values : List[float]
            List of x values.
        value_to_approximate : float
            The x value for which to perform interpolation.
        function : Union[None, Callable], optional
            The function to use for interpolation.
        function_values : Union[None, List[float]], optional
            List of function values corresponding to given_values.
        """
        super().__init__(given_values, value_to_approximate, function, function_values)

    def _solve(self, give_splines: bool = False):
        """
        Calculate linear spline solution.

        Parameters
        ----------
        give_splines : bool, optional
            Flag indicating whether to return spline values.

        Returns
        -------
            Tuple containing x values, y values, and spline values (if ``give_splines`` is True).
        """
        given, func_vals = self.given_values, self.function_values
        x_to_approx, data_points = self.value_to_approximate, len(self.given_values) - 1

        xs = [(given[i], given[i + 1]) for i in range(data_points)]
        ys = [(func_vals[i], func_vals[i + 1]) for i in range(data_points)]

        spline = [((x_to_approx - x2) * y1 - (x_to_approx - x1) * y2) / (x1 - x2) for (x1, x2), (y1, y2) in zip(xs, ys)]

        return (xs, ys, spline) if give_splines else spline

    def solution_set(self, give_splines: bool = False):
        """
        Get the solution set.

        Parameters
        ----------
        give_splines : bool, optional
            Flag indicating whether to return spline values.

        Returns
        -------
            Tuple containing x values, y values, and spline values (if ``give_splines`` is True).
        """
        return self._solve(give_splines=give_splines)

    def interpolate(self) -> float:
        """
        Perform linear spline interpolation.

        Returns
        -------
        float
            Interpolated value.
        """
        return get_solution(self.given_values, self.value_to_approximate, self._solve())

    def show_splines(self, full: bool = False):
        """
        Display linear spline equations.

        Parameters
        ----------
        full : bool, optional
            Flag indicating whether to display the full spline equations.
        """
        xs, ys, spline = self._solve(give_splines=True)

        print('The splines are approximated to 4 decimal places for display purposes only.\n')
        for i in range(len(spline)):
            den1 = xs[i][0] - xs[i][1]
            den2 = -den1
            if full:
                print(f'Sp{i + 1}: {ys[i][0] / den1:+}(x - {xs[i][1]}) {ys[i][1] / den2:+}(x - {xs[i][0]})')
            else:
                f1 = ys[i][0] / den1
                f1 += ys[i][1] / den2

                f2 = (ys[i][0] * xs[i][1]) / den1
                f2 += (ys[i][1] * xs[i][0]) / den2
                f2 *= -1

                print(f'Sp{i + 1}: {f1:+.4f}x {f2:+.4f}')
