"""Created on Nov 01 23:03:42 2023"""
import numpy as np

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


class QuadraticSpline(INTERPOLATION):

    def __init__(self, given_values, value_to_approximate, function=None, function_values=None, last_equation='first'):
        super().__init__(given_values, value_to_approximate, function, function_values)
        self.last_equation = last_equation

    @property
    def matrix_dimensions(self):
        n = len(self.given_values)
        return 3 * (n - 1)

    def __matrix_solver(self):
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

        return solution

    @property
    def solution_set(self):
        return self.__matrix_solver()

    def interpolate(self, n_derivative=0):
        approximation = self.value_to_approx
        solution = self.solution_set

        self.given_values.insert(0, self.value_to_approx)
        self.given_values.sort()
        idx_ = self.given_values.index(self.value_to_approx) - 1
        del self.given_values[idx_ + 1]

        req_solution = solution[idx_]

        if n_derivative == 0:
            return req_solution[0] * approximation**2 + req_solution[1] * approximation + req_solution[2]
        elif n_derivative == 1:
            return 2 * req_solution[0] * approximation + req_solution[1]
        elif n_derivative == 2:
            return 2 * req_solution[0]

    def show_splines(self) -> None:
        solution = self.solution_set

        print('The splines are approximated to 4 decimal places for display purposes only.\n')
        for i in range(len(solution)):
            print(f'Sp{i + 1}: {solution[i][0]:+.4f}x^2 {solution[i][1]:+.4f}x {solution[i][2]:+.4f} = 0')

        return None

# x = [0, 10, 15, 20, 22.5, 30]
# y = [0, 227.04, 362.78, 517.35, 602.97, 901.67]
# value_to_approx = 16
#
#
# def evaluate_spline(_x, solution_):
#     return solution_[0] * _x**2 + solution_[1] * _x + solution_[2]
#
#
# c1 = QuadraticSpline(x, value_to_approx, function_values=y, last_equation='first')
# ss = c1.solution_set
#
# x1 = np.linspace(x[0], x[1], 1000)
# x2 = np.linspace(x[1], x[2], 1000)
# x3 = np.linspace(x[2], x[3], 1000)
# x4 = np.linspace(x[3], x[4], 1000)
# x5 = np.linspace(x[4], x[5], 1000)
#
# plt.plot(x1, evaluate_spline(x1, ss[0]), 'b-')
# plt.plot(x2, evaluate_spline(x2, ss[1]), 'g-')
# plt.plot(x3, evaluate_spline(x3, ss[2]), 'r-')
# plt.plot(x4, evaluate_spline(x4, ss[3]), 'm-')
# plt.plot(x5, evaluate_spline(x5, ss[4]), 'y-')
#
# plt.plot(x, y, 'ko')
#
# plt.tight_layout()
# plt.show()
