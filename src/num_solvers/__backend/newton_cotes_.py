"""Newton-Cotes integration base module.

This module provides basic integration functionality of various Newton-Cotes integrators. Along with that, the module
provides a class for pretty-printing the integrator tables. The classes used in this module are,

- DataTable: Provides the pretty printing of integrator data tables.
- Iterators: An enumeration class for the integrators.
- SingleVariableNewtonCotesIntegrators: The base class for NC integrators.

Along with these classes, this module also provides a helper function to the newton_cotes_integrators front-end module

- solve_helper: Wraps the solver into a function with solution and data table.

Created on Nov 30 16:37:17 2023
"""

__all__ = ['SingleVariableNewtonCotesIntegrators', 'solve_helper']

import enum

from custom_inherit import doc_inherit

from .core_helpers_ import round_list_
from .. import DOC_STYLE, FList, Func, IFloat, N_DECIMAL


class DataTable:
    """Class for making integrator data tables."""

    def __init__(self, x: FList, y: FList):
        """
        Initializer for DataTable class

        Parameters
        ----------
        x:
            List of x values.
        y:
            List of y values.
        """
        self.x = x
        self.y = y

    def print_table(self):
        """Prints the data table."""

        # Find the maximum length of strings in x and y
        n_len = 6
        max_len_x = max(len(str(element)) for element in self.x)
        max_len_y = max(len(str(element)) for element in self.y)

        # Add extra space for padding
        max_len_x += 10
        max_len_y += 10
        n_len += 2

        print('-' * (n_len + max_len_x + max_len_y + 4))
        # Print the table header
        print(f"|{'n':^{n_len}}|{'x':^{max_len_x}}|{'y':^{max_len_y}}|")

        # Print the table border
        print('-' * (n_len + max_len_x + max_len_y + 4))

        # Print the table rows
        for index, value in enumerate(self.x):
            print(f"|{str(index + 1):^{n_len}}|{str(value):^{max_len_x}}|"
                  f"{str(self.y[index]):^{max_len_y}}|")

        # Print the table border
        print('-' * (n_len + max_len_x + max_len_y + 4))


@enum.unique
class Iterators(enum.Enum):
    """Iterator enum class."""

    TRAPEZOID = 1
    SIMPSON_13 = 2
    SIMPSON_38 = 3
    BOOLE = 4
    WEDDLE_3H10 = 5.1
    WEDDLE_41H140 = 5.2


class SingleVariableNewtonCotesIntegrators:
    """Provides functionality to integrate single variables using Newton-Cotes methods."""

    _I = Iterators

    str_iters = {'trapezoid': _I.TRAPEZOID,
                 'simpson_13': _I.SIMPSON_13,
                 'simpson_38': _I.SIMPSON_38,
                 'boole': _I.BOOLE,
                 'weddle_3h10': _I.WEDDLE_3H10,
                 'weddle_41h140': _I.WEDDLE_41H140}

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

    def __init__(self, function: Func, x_0: IFloat, x_n: IFloat, composites: int = 1,
                 integration_method: str = 'weddle'):
        """
        Initializer for Newton-Cotes single variable integrator class.

        Parameters
        ----------
        function :
            The function to be integrated.
        x_0 :
            The lower limit of integration.
        x_n :
            The upper limit of integration.
        composites :
            Number of sub-intervals for composite method. Default is 1.
        integration_method:
            The integration method to use for evaluation.
        """
        self.function = function
        self.x_0 = x_0
        self.x_n = x_n
        self.solver = self.str_iters[integration_method] if isinstance(integration_method, str) else integration_method

        _solver = self.solver

        self.composites = composites * self.nodes[_solver]
        self.dx = (x_n - x_0) / self.composites

        self.n_points = self.nodes[_solver]
        self.norm_value = self.norm[_solver](self.dx)
        self.outer_multiplier_value = self.outer_multipliers[_solver]
        self.inner_multipliers_values = self.inner_multipliers[_solver]

    def _value_table(self):
        x, dx, func, composites = self.x_0, self.dx, self.function, self.composites

        x_values = [x + i * dx for i in range(composites + 1)]
        y_values = [func(i) for i in x_values]

        return x_values, y_values

    def value_table(self, n_decimal: int = N_DECIMAL):
        """
        Prints the data table for the given values.

        Parameters
        ----------
        n_decimal:
            The number of decimal places to round up to. Default is 8.
        """

        vals = self._value_table()

        x_values = round_list_(vals[0], n_decimal)
        y_values = round_list_(vals[1], n_decimal)

        DataTable(x_values, y_values).print_table()

    def solve(self) -> IFloat:
        """
        Solves the Newton-Cotes integrator.

        Returns
        -------
            The integration value for the provided function.
        """

        value_table, solver_method, normalizer = self._value_table()[1], self.solver, self.norm_value
        inner_mul, outer_mul, n_points = self.inner_multipliers_values, self.outer_multiplier_value, self.n_points

        f_0 = [outer_mul * i for i in [value_table[0], value_table[-1]]]

        f_list = []
        for multiplier, start in zip(inner_mul, range(1, n_points + 1)):
            if solver_method == self._I.TRAPEZOID:
                f_list.extend(multiplier * i for i in value_table[start:-1])
            else:
                f_list.extend(multiplier * i for i in value_table[start:-1:n_points])

        f_list[0:0] = f_0

        return normalizer * sum(f_list)


@doc_inherit(SingleVariableNewtonCotesIntegrators.__init__, style=DOC_STYLE)
def solve_helper(function: Func, x_0: IFloat, x_n: IFloat, composites: int = 1, get_table: bool = False,
                 method: str = 'trapezoid'):
    """
    Helper function for performing numerical integration.

    Returns
    -------
        The numerical result of integration.
    """

    solver = SingleVariableNewtonCotesIntegrators(function, x_0, x_n, composites, method)

    if get_table:
        solver.value_table()

    return solver.solve()
