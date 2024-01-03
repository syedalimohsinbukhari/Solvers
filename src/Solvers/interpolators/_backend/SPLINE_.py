"""Created on Dec 28 17:28:17 2023"""

from typing import List

import numpy as np

from . import TOLERANCE
from .INTERPOLATION_ import INTERPOLATION


# TODO: Add docstring


class SPLINE(INTERPOLATION):

    def __init__(self, given_values, value_to_approximate, function=None, function_values=None):
        super().__init__(given_values, value_to_approximate, function, function_values)

    def _solve(self):
        pass

    def show_splines(self):
        pass

    def solution_set(self):
        pass


class LinearSpline(SPLINE):

    def __init__(self, given_values, value_to_approximate, function=None, function_values=None) -> None:
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

    def _solve(self, give_values: bool = False):
        """
        Calculate linear spline solution.

        Parameters
        ----------
        give_values : bool, optional
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

        return (xs, ys, spline) if give_values else spline

    def solution_set(self, give_values: bool = False):
        """
        Get the solution set.

        Parameters
        ----------
        give_values : bool, optional
            Flag indicating whether to return spline values.

        Returns
        -------
            Tuple containing x values, y values, and spline values (if ``give_splines`` is True).
        """
        return self._solve(give_values=give_values)

    def interpolate(self) -> float:
        """
        Perform linear spline interpolation.

        Returns
        -------
        float
            Interpolated value.
        """
        return get_solution(self.given_values, self.value_to_approximate, self._solve())

    def show_splines(self):
        """
        Display linear spline equations.

        """
        xs, ys, spline = self._solve(give_values=True)

        print('The splines are approximated to 4 decimal places for display purposes only.')
        for i in range(len(spline)):
            den1 = xs[i][0] - xs[i][1]
            den2 = -den1

            f1 = ys[i][0] / den1
            f1 += ys[i][1] / den2

            f2 = (ys[i][0] * xs[i][1]) / den1
            f2 += (ys[i][1] * xs[i][0]) / den2
            f2 *= -1

            print(f'Sp{i + 1}: {f1:+08.4f}x {f2:+08.4f}'
                  f'\t; x ∈ {"[" if i == 0 else "("}{self.given_values[i]:+}, {self.given_values[i + 1]:+}]')


class QuadraticSpline(SPLINE):

    def __init__(self, given_values, value_to_approximate, function=None, function_values=None, last_equation='first'):
        super().__init__(given_values, value_to_approximate, function, function_values)
        self.last_equation = last_equation

    @property
    def matrix_dimensions(self):
        n = len(self.given_values)
        return 3 * (n - 1)

    def _solve(self):
        len_mat = self.matrix_dimensions
        n_given = len(self.given_values)

        matrix_a = [[0] * len_mat for _ in range(len_mat)]
        matrix_b = [0] * len_mat

        i_index = spline_equations(self.given_values, self.function_values, matrix_a, matrix_b, 0, 2)
        combined = middle_point_equations(self.given_values, 2)
        fill_middle_point_splines(combined, i_index, matrix_a, n_given, 2)

        if self.last_equation == 'first':
            matrix_a[-1][0] = 1
        elif self.last_equation == 'last':
            matrix_a[-1][-3] = 1
        else:
            raise ValueError('Possible values are \'first\' or \'last\'.')

        solution = apply_linalg(matrix_a, matrix_b, n_given, 2)

        return solution

    def interpolate(self, n_derivative=0):
        approximation = self.value_to_approximate
        req_solution = get_solution(self.given_values, self.value_to_approximate, self._solve())

        if n_derivative == 0:
            return req_solution[0] * approximation**2 + req_solution[1] * approximation + req_solution[2]
        elif n_derivative == 1:
            return 2 * req_solution[0] * approximation + req_solution[1]
        elif n_derivative == 2:
            return 2 * req_solution[0]
        else:
            return 0

    def solution_set(self, get_matrix=False):
        return self._solve()

    def show_splines(self):
        solution = self._solve()

        print('The splines are approximated to 4 decimal places for display purposes only.')
        for i in range(len(solution)):
            print(f'Sp{i + 1}: {solution[i][0]:+08.4f}x^2 {solution[i][1]:+08.4f}x {solution[i][2]:+08.4f}'
                  f'\t; x ∈ {"[" if i == 0 else "("}{self.given_values[i]:+}, {self.given_values[i + 1]:+}]')


class NaturalCubicSpline(SPLINE):

    def __init__(self, given_values, value_to_approximate, function=None, function_values=None, last_equation='first'):
        super().__init__(given_values, value_to_approximate, function, function_values)
        self.last_equation = last_equation

    @property
    def matrix_dimensions(self):
        n = len(self.given_values)
        return 4 * (n - 1)

    def _solve(self):
        len_mat = self.matrix_dimensions
        n_given = len(self.given_values)

        matrix_a = [[0] * len_mat for _ in range(len_mat)]
        matrix_b = [0] * len_mat

        i_index1 = spline_equations(self.given_values, self.function_values, matrix_a, matrix_b, 0, 3)
        combined = middle_point_equations(self.given_values, 3)
        combined_der = middle_point_equations(self.given_values, -1)

        i_index1 = fill_middle_point_splines(combined, i_index1, matrix_a, n_given, 3)
        fill_middle_point_splines(combined_der, i_index1, matrix_a, n_given, 3)

        matrix_a[-2][0] = 6 * matrix_a[0][2]
        matrix_a[-2][1] = 2 * matrix_a[0][3]
        matrix_a[-1][-4] = 6 * matrix_a[2 * (len(self.given_values) - 1) - 1][-2]
        matrix_a[-1][-3] = 2 * matrix_a[2 * (len(self.given_values) - 1) - 1][-1]

        solution = apply_linalg(matrix_a, matrix_b, n_given, 3)

        return solution

    def solution_set(self, get_matrix=False):
        return self._solve()

    def show_splines(self):
        solution = self._solve()

        print('The splines are approximated to 4 decimal places for display purposes only.')
        for i in range(len(solution)):
            print(f'Sp{i + 1}: {solution[i][0]:+08.4f}x^3 {solution[i][1]:+08.4f}x^2 {solution[i][2]:+08.4f}x '
                  f'{solution[i][3]:+08.4f}\t; x ∈ {"[" if i == 0 else "("}{self.given_values[i]:+}, '
                  f'{self.given_values[i + 1]:+}]')

    def interpolate(self, n_derivative=0):
        approx = self.value_to_approximate
        solution = get_solution(self.given_values, self.value_to_approximate, self._solve())

        if n_derivative == 0:
            return solution[0] * approx**3 + solution[1] * approx**2 + solution[2] * approx + solution[3]
        elif n_derivative == 1:
            return 3 * solution[0] * approx**2 + 2 * solution[1] * approx + solution[2]
        elif n_derivative == 2:
            return 6 * solution[0] * approx + 2 * solution[1]
        elif n_derivative == 3:
            return 6 * solution[0]
        else:
            return 0


def fill_middle_point_splines(combined, i_index1, matrix_a, n_given, degree_of_spline):
    degree_of_spline, j_index, k_index = degree_of_spline + 1, 0, 0
    i_temp = i_index1
    while i_index1 < i_temp + n_given - 2:
        ind = degree_of_spline * k_index
        matrix_a[i_index1][ind:ind + degree_of_spline] = combined[0][j_index]
        matrix_a[i_index1][ind + degree_of_spline:ind + 2 * degree_of_spline] = combined[1][j_index]

        i_index1 += 1
        j_index += 1
        k_index += 1
    return i_index1


def spline_equations(given_values, function_values, matrix_a, matrix_b, i_index, degree_of_spline):
    cond, degree_of_spline = 2 * (len(given_values) - 1), degree_of_spline + 1
    j_index, k_index = 0, 0

    while i_index < cond:
        ind = degree_of_spline * k_index
        val1, val2 = given_values[j_index], given_values[j_index + 1]

        if degree_of_spline == 3:
            matrix_a[i_index][ind:ind + degree_of_spline] = [val1**2, val1, 1]
            matrix_a[i_index + 1][ind:ind + degree_of_spline] = [val2**2, val2, 1]
        else:
            matrix_a[i_index][ind:ind + degree_of_spline] = [val1**3, val1**2, val1, 1]
            matrix_a[i_index + 1][ind:ind + degree_of_spline] = [val2**3, val2**2, val2, 1]

        matrix_b[i_index] = function_values[j_index]
        matrix_b[i_index + 1] = function_values[j_index + 1]

        i_index += 2
        j_index += 1
        k_index += 1

    return i_index


def middle_point_equations(given_values, degree_of_spline):
    values, degree_of_spline = given_values[1:-1], degree_of_spline + 1

    if degree_of_spline == 3:
        a_ = [[2 * j * i for i in values] for j in [1, -1]]
        b_ = [[j for _ in values] for j in [1, -1]]
        c_ = [[0] * len(values) for _ in range(2)]

        res = [a_, b_, c_]
    elif degree_of_spline == 4:
        a_ = [[3 * j * i**2 for i in values] for j in [1, -1]]
        b_ = [[2 * j * _ for _ in values] for j in [1, -1]]
        c_ = [[j for _ in values] for j in [1, -1]]
        d_ = [[0] * len(values) for _ in range(2)]

        res = [a_, b_, c_, d_]
    else:
        a_ = [[6 * j * i for i in values] for j in [1, -1]]
        b_ = [[2 * j for _ in values] for j in [1, -1]]
        c_ = [[0] * len(values) for _ in range(2)]
        d_ = c_

        res = [a_, b_, c_, d_]

    return [[[*p] for p in zip(*_)] for _ in zip(*res)]


def apply_linalg(matrix_a, matrix_b, n_given, degree_of_spline):
    # solve using numpy solver
    solution = np.linalg.solve(np.array(matrix_a), np.array(matrix_b))
    # reshape according to the data
    solution = np.reshape(solution, (n_given - 1, degree_of_spline + 1))
    # convert all very small values to 0
    solution[np.abs(solution) < TOLERANCE] = 0

    return solution


def get_solution(given_values, value_to_approximate, solution):
    given_values.insert(0, value_to_approximate)
    given_values.sort()
    idx_ = given_values.index(value_to_approximate) - 1
    del given_values[idx_ + 1]

    return solution[idx_]
