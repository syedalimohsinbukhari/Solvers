"""Created on Nov 01 23:03:42 2023"""
import numpy as np
from matplotlib import pyplot as plt

from src.interpolators.interpolation.INTERPOLATION_ import INTERPOLATION


def get_solution(given_values, value_to_approximate, solution):
    given_values.insert(0, value_to_approximate)
    given_values.sort()
    idx_ = given_values.index(value_to_approx) - 1
    del given_values[idx_ + 1]

    return solution[idx_]


class LinearSpline(INTERPOLATION):

    def __init__(self, given_values, value_to_approximate, function=None, function_values=None):
        super().__init__(given_values, value_to_approximate, function, function_values)

    def __spline_solver(self, give_splines=False):
        x_to_approx = self.value_to_approx
        data_points = len(self.given_values)

        xs = [[self.given_values[i], self.given_values[i + 1]] for i in range(data_points - 1)]
        ys = [[self.function_values[i], self.function_values[i + 1]] for i in range(data_points - 1)]

        spline = [((x_to_approx - x2) * y1 - (x_to_approx - x1) * y2) / (x1 - x2) for (x1, x2), (y1, y2) in zip(xs, ys)]

        return xs, ys, spline if give_splines else spline

    def solution_set(self, give_splines=False):
        return self.__spline_solver(give_splines=give_splines)

    def interpolate(self):
        return get_solution(self.given_values, self.value_to_approx, self.solution_set())

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

    def plot_splines(self):
        given = self.given_values
        xs, ys, _ = self.solution_set(give_splines=True)

        def get_spline(x_par, y_par, point):
            num1, num2 = (point - x_par[1]) * y_par[0], (point - x_par[0]) * y_par[1]
            denominator = x_par[0] - x_par[1]

            return (num1 - num2) / denominator

        def get_xs(st, ed, n=10):
            return np.linspace(st, ed, n)

        x_vals = [get_xs(i, j, 100) for i, j in zip(given[1:], given[:-1])]
        y_vals = [get_spline(xs[i], ys[i], v) for i, v in enumerate(x_vals)]

        plt.figure()
        plt.plot(given, self.function_values, 'k--')
        [plt.plot(i, j, ls='--') for i, j in zip(x_vals, y_vals)]
        plt.plot(given, self.function_values, 'ko')
        plt.title('Linear Spline Approximation.')
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.tight_layout()
        plt.show()


class QuadraticSpline(INTERPOLATION):

    def __init__(self, given_values, value_to_approximate, function=None, function_values=None, last_equation='first'):
        super().__init__(given_values, value_to_approximate, function, function_values)
        self.last_equation = last_equation

    @property
    def matrix_dimensions(self):
        n = len(self.given_values)
        return 3 * (n - 1)

    def __matrix_solver(self, give_matrix=False):
        len_mat = self.matrix_dimensions
        n_given = len(self.given_values)
        matrix_a = [[0 for _ in range(len_mat)] for _ in range(len_mat)]
        matrix_b = [0 for _ in range(len_mat)]

        i, j, k = 0, 0, 0
        while i < len_mat - n_given + 1:
            ind = 3 * k
            val1, val2 = self.given_values[j], self.given_values[j + 1]
            matrix_a[i][ind:ind + 3] = [val1**2, val1, 1]
            matrix_a[i + 1][ind:ind + 3] = [val2**2, val2, 1]

            matrix_b[i] = self.function_values[j]
            matrix_b[i + 1] = self.function_values[j + 1]

            i += 2
            j += 1
            k += 1

        given_values = self.given_values[1:-1]
        a_ = [[2 * j * i for i in given_values] for j in [1, -1]]
        b_ = [[j for _ in given_values] for j in [1, -1]]
        c_ = [[0] * len(given_values) for _ in range(2)]

        combined = [[[*p] for p in zip(x1, y1, z1)] for x1, y1, z1 in zip(a_, b_, c_)]

        j, k = 0, 0
        while i < len_mat - 1:
            ind = 3 * k
            matrix_a[i][ind:ind + 3] = combined[0][j]
            matrix_a[i][ind + 3:ind + 3 + 3] = combined[1][j]

            i += 1
            j += 1
            k += 1

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
        return self.__matrix_solver(give_matrix=give_matrix)

    def interpolate(self, n_derivative=0):
        approximation = self.value_to_approx
        req_solution = get_solution(self.given_values, self.value_to_approx, self.solution_set())

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

    def plot_splines(self):

        def get_splines(pars, _x):
            return pars[0] * _x**2 + pars[1] * _x + pars[2]

        def get_xs(st, ed, n=10):
            return np.linspace(st, ed, n)

        vals = self.given_values
        ss = self.solution_set()

        xs = [get_xs(i, j) for i, j in zip(vals[:-1], vals[1:])]
        splines = [get_splines(i, j) for i, j in zip(ss, xs)]

        plt.figure()
        plt.plot(vals, self.function_values, 'k--')
        [plt.plot(i, j) for i, j in zip(xs, splines)]
        plt.plot(vals, self.function_values, 'ko')
        plt.title('Quadratic Spline Approximation.')
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.tight_layout()
        plt.show()


def f_of_x(_x):
    return _x**2 - (np.tan(_x**5) / _x)**-1


# x1 = [2.2, 8.2, 12.2, 20.2, 21.5, 28.9, 35, 45]
# y1 = [f_of_x(i) for i in x1]

x1 = [0, 10, 15, 20, 22.5, 30]
y1 = [0, 227.04, 362.78, 517.35, 602.97, 901.67]
value_to_approx = 16
c1 = LinearSpline(x1, value_to_approx, function_values=y1)
ssc1 = c1.solution_set()

c2 = QuadraticSpline(x1, value_to_approx, function_values=y1, last_equation='last')
ssc2 = c2.solution_set()


def resolve_x(start, end, num=100):
    return np.linspace(start, end, num)


def quad_spline(pars, _x):
    return pars[0] * _x**2 + pars[1] * _x + pars[2]


def linear_spline(x_par, y_par, point):
    num1, num2 = (point - x_par[1]) * y_par[0], (point - x_par[0]) * y_par[1]
    denominator = x_par[0] - x_par[1]

    return (num1 - num2) / denominator


# x_resolved = [resolve_x(i, j, 100) for i, j in zip(x1[1:], x1[:-1])]
# lin_ = [linear_spline(ssc1[0][i], ssc1[1][i], v) for i, v in enumerate(x_resolved)]
# quad_ = [quad_spline(i, j) for i, j in zip(ssc2, x_resolved)]
#
# plt.plot(x1, y1, 'k--')
# [plt.plot(i, j, color='b') for i, j in zip(x_resolved, lin_)]
# [plt.plot(i, j, color='r') for i, j in zip(x_resolved, quad_)]
# plt.plot(x1, y1, 'ko')
# plt.show()

c1.show_splines()
c1.plot_splines()
