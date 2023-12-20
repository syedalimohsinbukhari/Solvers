"""Created on Nov 01 23:03:42 2023"""

import numpy as np

from .interpolation.INTERPOLATION_ import INTERPOLATION


# TODO: Clean up the spline code and see if it can be generalized

class SPLINE(INTERPOLATION):

    def __init__(self, given_values, value_to_approximate, function=None, function_values=None):
        super().__init__(given_values, value_to_approximate, function, function_values)

    def _solve(self):
        pass

    def show_splines(self):
        pass

    @staticmethod
    def x_values(start, end, n=1000):
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

    def __init__(self, given_values, value_to_approximate, function=None, function_values=None):
        super().__init__(given_values, value_to_approximate, function, function_values)

    def _solve(self, give_splines=False):
        x_to_approx = self.value_to_approximate
        data_points = len(self.given_values)

        xs = [[self.given_values[i], self.given_values[i + 1]] for i in range(data_points - 1)]
        ys = [[self.function_values[i], self.function_values[i + 1]] for i in range(data_points - 1)]

        spline = [((x_to_approx - x2) * y1 - (x_to_approx - x1) * y2) / (x1 - x2) for (x1, x2), (y1, y2) in zip(xs, ys)]

        return xs, ys, spline if give_splines else spline

    def solution_set(self, give_splines=False):
        return self._solve(give_splines=give_splines)

    def interpolate(self):
        return get_solution(self.given_values, self.value_to_approximate, self.solution_set())

    def show_splines(self, full=False) -> None:
        xs, ys, spline = self.solution_set(give_splines=True)

        print('The splines are approximated to 4 decimal places for display purposes only.\n')
        for i in range(len(spline)):
            den1 = xs[i][0] - xs[i][1]
            den2 = -den1
            if full:
                print(f'Sp{i + 1}: {ys[i][0] / den1}(x - {xs[i][1]}) - {ys[i][1] / den2}(x - {xs[i][0]})')
            else:
                print(f'Sp{i + 1}: {ys[i][0] / den1:+.4f}(x - {xs[i][1]}) {ys[i][1] / den2:+.4f}(x - {xs[i][0]})')

        return None


class QuadraticSpline(SPLINE):

    def __init__(self, given_values, value_to_approximate, function=None, function_values=None, last_equation='first'):
        super().__init__(given_values, value_to_approximate, function, function_values)
        self.last_equation = last_equation

    @property
    def matrix_dimensions(self):
        n = len(self.given_values)
        return 3 * (n - 1)

    def _solve(self, give_matrix=False):
        len_mat = self.matrix_dimensions
        n_given = len(self.given_values)
        matrix_a = [[0 for _ in range(len_mat)] for _ in range(len_mat)]
        matrix_b = [0 for _ in range(len_mat)]

        i_index, j_index, k_index = 0, 0, 0
        cond = 2 * (n_given - 1)
        while i_index < cond:
            ind = 3 * k_index
            val1, val2 = self.given_values[j_index], self.given_values[j_index + 1]
            matrix_a[i_index][ind:ind + 3] = [val1**2, val1, 1]
            matrix_a[i_index + 1][ind:ind + 3] = [val2**2, val2, 1]

            matrix_b[i_index] = self.function_values[j_index]
            matrix_b[i_index + 1] = self.function_values[j_index + 1]

            i_index += 2
            j_index += 1
            k_index += 1

        given_values = self.given_values[1:-1]
        a_ = [[2 * j * i for i in given_values] for j in [1, -1]]
        b_ = [[j for _ in given_values] for j in [1, -1]]
        c_ = [[0] * len(given_values) for _ in range(2)]

        combined = [[[*p] for p in zip(x1, y1, z1)] for x1, y1, z1 in zip(a_, b_, c_)]

        j_index, k_index = 0, 0
        while i_index < len_mat - 1:
            ind = 3 * k_index
            matrix_a[i_index][ind:ind + 3] = combined[0][j_index]
            matrix_a[i_index][ind + 3:ind + 3 + 3] = combined[1][j_index]

            i_index += 1
            j_index += 1
            k_index += 1

        if self.last_equation == 'first':
            matrix_a[-1][0] = 1
        elif self.last_equation == 'last':
            matrix_a[-1][-3] = 1
        else:
            raise ValueError('Possible values are \'first\' or \'last\'.')

        # solve using numpy solver
        solution = np.linalg.solve(np.array(matrix_a), np.array(matrix_b))
        # reshape according to the data
        solution = np.reshape(solution, (n_given - 1, 3))
        # convert all very small values to 0
        solution[np.abs(solution) < 1e-10] = 0

        return [matrix_a, solution] if give_matrix else solution

    def solution_set(self, give_matrix=False):
        return self._solve(give_matrix=give_matrix)

    def interpolate(self, n_derivative=0):
        approximation = self.value_to_approximate
        req_solution = get_solution(self.given_values, self.value_to_approximate, self.solution_set())

        if n_derivative == 0:
            return req_solution[0] * approximation**2 + req_solution[1] * approximation + req_solution[2]
        elif n_derivative == 1:
            return 2 * req_solution[0] * approximation + req_solution[1]
        elif n_derivative == 2:
            return 2 * req_solution[0]

    def show_splines(self, full=False) -> None:
        solution = self.solution_set()

        print('The splines are approximated to 4 decimal places for display purposes only.\n')
        for i in range(len(solution)):
            if full:
                print(f'Sp{i + 1}: {solution[i][0]}x^2 {solution[i][1]}x {solution[i][2]} = 0')
            else:
                print(f'Sp{i + 1}: {solution[i][0]:+.4f}x^2 {solution[i][1]:+.4f}x {solution[i][2]:+.4f} = 0')

        return None


class NaturalCubicSpline(SPLINE):

    def __init__(self, given_values, value_to_approximate, function=None, function_values=None, last_equation='first'):
        super().__init__(given_values, value_to_approximate, function, function_values)
        self.last_equation = last_equation

    @property
    def matrix_dimensions(self):
        n = len(self.given_values)
        return 4 * (n - 1)

    def _solve(self, give_matrix=False):
        len_mat = self.matrix_dimensions
        n_given = len(self.given_values)
        matrix_a = [[0 for _ in range(len_mat)] for _ in range(len_mat)]
        matrix_b = [0 for _ in range(len_mat)]

        cond = 2 * (n_given - 1)
        i_index1, j_index, k_index = 0, 0, 0
        while i_index1 < cond:
            ind = 4 * k_index
            val1, val2 = self.given_values[j_index], self.given_values[j_index + 1]
            matrix_a[i_index1][ind:ind + 4] = [val1**3, val1**2, val1, 1]
            matrix_a[i_index1 + 1][ind:ind + 4] = [val2**3, val2**2, val2, 1]

            matrix_b[i_index1] = self.function_values[j_index]
            matrix_b[i_index1 + 1] = self.function_values[j_index + 1]

            i_index1 += 2
            j_index += 1
            k_index += 1

        given_values = self.given_values[1:-1]
        a_ = [[3 * j * i**2 for i in given_values] for j in [1, -1]]
        b_ = [[2 * j * _ for _ in given_values] for j in [1, -1]]
        c_ = [[j for _ in given_values] for j in [1, -1]]
        d_ = [[0] * len(given_values) for _ in range(2)]
        combined = [[[*p] for p in zip(w1, x1, y1, z1)] for w1, x1, y1, z1 in zip(a_, b_, c_, d_)]

        a_der = [[6 * j * i for i in given_values] for j in [1, -1]]
        b_der = [[2 * j for _ in given_values] for j in [1, -1]]
        c_der = [[0] * len(given_values) for _ in range(2)]
        d_der = c_der
        combined_der = [[[*p] for p in zip(w2, x2, y2, z2)] for w2, x2, y2, z2 in zip(a_der, b_der, c_der, d_der)]

        j_index, k_index = 0, 0
        i_temp = i_index1
        while i_index1 < i_temp + n_given - 2:
            ind = 4 * k_index
            matrix_a[i_index1][ind:ind + 4] = combined[0][j_index]
            matrix_a[i_index1][ind + 4:ind + 4 + 4] = combined[1][j_index]

            i_index1 += 1
            j_index += 1
            k_index += 1

        j_index, k_index = 0, 0
        i_temp = i_index1
        while i_index1 < i_temp + n_given - 2:
            ind = 4 * k_index
            matrix_a[i_index1][ind:ind + 4] = combined_der[0][j_index]
            matrix_a[i_index1][ind + 4:ind + 4 + 4] = combined_der[1][j_index]

            i_index1 += 1
            j_index += 1
            k_index += 1

        matrix_a[-2][0] = 6
        matrix_a[-2][1] = 1
        matrix_a[-1][-4] = 6
        matrix_a[-1][-3] = 1

        # solve using numpy solver
        solution = np.linalg.solve(np.array(matrix_a), np.array(matrix_b))
        # reshape according to the data
        solution = np.reshape(solution, (n_given - 1, 4))
        # convert all very small values to 0
        solution[np.abs(solution) < 1e-10] = 0

        return [matrix_a, solution] if give_matrix else solution

    def solution_set(self, give_matrix=False):
        return self._solve(give_matrix=give_matrix)

    def show_splines(self, full=False) -> None:
        solution = self.solution_set()

        print('The splines are approximated to 4 decimal places for display purposes only.\n')
        for i in range(len(solution)):
            if full:
                print(f'Sp{i + 1}: {solution[i][0]}x^3 {solution[i][1]}x^2 {solution[i][2]}x + {solution[i][3]} = 0')
            else:
                print(f'Sp{i + 1}: {solution[i][0]:+.4f}x^3 {solution[i][1]:+.4f}x^2 {solution[i][2]:+.4f}x '
                      f'{solution[i][3]:+.4f} = 0')

        return None
