"""Created on Nov 29 14:48:53 2023"""

from enum import Enum


def newton_raphson_solver(function, func_derivative, x_0, tol=1e-14):
    f, df = function, func_derivative

    if abs(f(x_0)) < tol:
        return x_0
    else:
        return newton_raphson_solver(f, df, x_0 - f(x_0) / df(x_0), tol)


class _DataTable:
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def print_table(self, precision=100):
        # Find the maximum length of strings in x and y
        n_len = 6
        max_len_x = max(len(str(element)) for element in self.x)
        max_len_y = max(len(str(element)) for element in self.y)

        # Add extra space for padding
        max_len_x += 4
        max_len_y += 4
        n_len += 2

        print('-' * (n_len + max_len_x + max_len_y + 4))
        # Print the table header
        print(f"|{'n':^{n_len}}|{'x':^{max_len_x}}|{'y':^{max_len_y}}|")

        # Print the table border
        print('-' * (n_len + max_len_x + max_len_y + 4))

        # Print the table rows
        for i in range(len(self.x)):
            print(f"|{str(i + 1):^{n_len}}|{str(self.x[i]):^{max_len_x}}|"
                  f"{str(round(self.y[i], precision)):^{max_len_y}}|")

        # Print the table border
        print('-' * (n_len + max_len_x + max_len_y + 4))


class Iterators(Enum):
    TRAPEZOID = 1
    SIMPSON_13 = 2
    SIMPSON_38 = 3
    BOOLE = 4
    WEDDLE_3H10 = 5.1
    WEDDLE_41H140 = 5.2


class NewtonCotesIterators:

    str_iters = {'trapezoid': Iterators.TRAPEZOID,
                 'simpson_13': Iterators.SIMPSON_13,
                 'simpson_38': Iterators.SIMPSON_38,
                 'boole': Iterators.BOOLE,
                 'weddle': Iterators.WEDDLE_3H10,
                 'weddle_long': Iterators.WEDDLE_41H140}

    nodes = {Iterators.TRAPEZOID: 2,
             Iterators.SIMPSON_13: 2,
             Iterators.SIMPSON_38: 3,
             Iterators.BOOLE: 4,
             Iterators.WEDDLE_3H10: 6,
             Iterators.WEDDLE_41H140: 6}

    norm = {Iterators.TRAPEZOID: lambda x: x / 2,
            Iterators.SIMPSON_13: lambda x: x / 3,
            Iterators.SIMPSON_38: lambda x: (3 * x) / 8,
            Iterators.BOOLE: lambda x: (2 * x) / 45,
            Iterators.WEDDLE_3H10: lambda x: (3 * x) / 10,
            Iterators.WEDDLE_41H140: lambda x: x / 140}

    outer_multipliers = {Iterators.TRAPEZOID: 1,
                         Iterators.SIMPSON_13: 1,
                         Iterators.SIMPSON_38: 1,
                         Iterators.BOOLE: 7,
                         Iterators.WEDDLE_3H10: 1,
                         Iterators.WEDDLE_41H140: 41}

    inner_multipliers = {Iterators.TRAPEZOID: [2],
                         Iterators.SIMPSON_13: [4, 2],
                         Iterators.SIMPSON_38: [3, 3, 2],
                         Iterators.BOOLE: [32, 12, 32, 2 * 7],
                         Iterators.WEDDLE_3H10: [5, 1, 6, 1, 5, 2 * 1],
                         Iterators.WEDDLE_41H140: [216, 27, 272, 27, 216, 2 * 41]}

    error_formula = {Iterators.TRAPEZOID: lambda x: x**2 / 12,
                     Iterators.SIMPSON_13: lambda x: x**4 / 180}

    def __init__(self, function, x_0, x_n, composites=1, solver='weddle'):
        self.function = function
        self.x_0 = x_0
        self.x_n = x_n
        self.solver = self.str_iters[solver]

        self.composites = composites * self.nodes[self.solver]
        self.dx = (x_n - x_0) / self.composites

    def _value_table(self):
        x_values = [self.x_0 + i * self.dx for i in range(self.composites + 1)]

        return {'x': x_values, 'y': [self.function(i) for i in x_values]}

    def value_table(self, precision=100):
        vals = self._value_table()
        _DataTable(vals['x'], vals['y']).print_table(precision=precision)

    def solve(self):
        n_points = self.nodes[self.solver]
        norm = self.norm[self.solver](self.dx)

        fs = self._value_table()['y']
        f_0 = [self.outer_multipliers[self.solver] * i for i in [fs[0], fs[-1]]]

        f_list = []
        for multiplier, start in zip(self.inner_multipliers[self.solver], range(1, n_points + 1)):
            if self.solver == Iterators.TRAPEZOID:
                f_list.extend(multiplier * i for i in fs[start:-1])
            else:
                f_list.extend(multiplier * i for i in fs[start:-1:n_points])

        f_list[0:0] = reversed(f_0)

        return norm * sum(f_list)
