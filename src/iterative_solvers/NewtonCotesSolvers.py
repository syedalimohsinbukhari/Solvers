"""Created on Nov 29 14:48:53 2023"""

from src.iterative_solvers.gaussian.NEWTON_COTES import Iterators, _DataTable


def newton_raphson_solver(function, func_derivative, x_0, tol=1e-14):
    f, df = function, func_derivative

    if abs(f(x_0)) < tol:
        return x_0
    else:
        return newton_raphson_solver(f, df, x_0 - f(x_0) / df(x_0), tol)


class NewtonCotesIterators:

    _I = Iterators

    str_iters = {'trapezoid': _I.TRAPEZOID,
                 'simpson_13': _I.SIMPSON_13,
                 'simpson_38': _I.SIMPSON_38,
                 'boole': _I.BOOLE,
                 'weddle': _I.WEDDLE_3H10,
                 'weddle_long': _I.WEDDLE_41H140}

    nodes = {_I.TRAPEZOID: 2, _I.SIMPSON_13: 2, _I.SIMPSON_38: 3, _I.BOOLE: 4, _I.WEDDLE_3H10: 6, _I.WEDDLE_41H140: 6}

    norm = {_I.TRAPEZOID: lambda x: x / 2,
            _I.SIMPSON_13: lambda x: x / 3,
            _I.SIMPSON_38: lambda x: (3 * x) / 8,
            _I.BOOLE: lambda x: (2 * x) / 45,
            _I.WEDDLE_3H10: lambda x: (3 * x) / 10,
            _I.WEDDLE_41H140: lambda x: x / 140}

    outer_multipliers = {_I.TRAPEZOID: 1, _I.SIMPSON_13: 1, _I.SIMPSON_38: 1, _I.BOOLE: 7, _I.WEDDLE_3H10: 1,
                         _I.WEDDLE_41H140: 41}

    inner_multipliers = {_I.TRAPEZOID: [2],
                         _I.SIMPSON_13: [4, 2],
                         _I.SIMPSON_38: [3, 3, 2],
                         _I.BOOLE: [32, 12, 32, 2 * 7],
                         _I.WEDDLE_3H10: [5, 1, 6, 1, 5, 2 * 1],
                         _I.WEDDLE_41H140: [216, 27, 272, 27, 216, 2 * 41]}

    def __init__(self, function, x_0, x_n, composites=1, solver='weddle'):
        self.function = function
        self.x_0 = x_0
        self.x_n = x_n
        self.solver = self.str_iters[solver] if isinstance(solver, str) else solver

        self.composites = composites * self.nodes[self.solver]
        self.dx = (x_n - x_0) / self.composites

        self.n_points = self.nodes[self.solver]
        self.norm_value = self.norm[self.solver](self.dx)
        self.outer_multiplier_value = self.outer_multipliers[self.solver]
        self.inner_multipliers_values = self.inner_multipliers[self.solver]

    def _value_table(self):
        x_values = [self.x_0 + i * self.dx for i in range(self.composites + 1)]
        return {'x': x_values, 'y': [self.function(i) for i in x_values]}

    def value_table(self, precision=100):
        vals = self._value_table()
        _DataTable(vals['x'], vals['y']).print_table(precision=precision)

    def solve(self):
        fs = self._value_table()['y']
        f_0 = [self.outer_multiplier_value * i for i in [fs[0], fs[-1]]]

        f_list = []
        for multiplier, start in zip(self.inner_multipliers_values, range(1, self.n_points + 1)):
            if self.solver == self._I.TRAPEZOID:
                f_list.extend(multiplier * i for i in fs[start:-1])
            else:
                f_list.extend(multiplier * i for i in fs[start:-1:self.n_points])

        f_list[0:0] = reversed(f_0)

        return self.norm_value * sum(f_list)
