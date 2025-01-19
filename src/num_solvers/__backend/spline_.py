"""Spline interpolation base module

This module provides the basic classes to perform different types of spline interpolations,
- SPLINE: The base class for spline interpolation.
- LinearSpline: For performing linear spline interpolation.
- QuadraticSpline: For performing quadratic spline interpolation.
- NaturalCubicSpline: For performing natural cubic spline interpolation.

The module also uses a number of private functions.
- _fill_initial_splines: To fill initial spline values in the matrix.
- _generate_middle_point_equations: To generates the equations for the middle points in the interpolation.
- _fill_middle_splines: To fill the generated middle point splines to the matrix.
- _matrix_solver: To solve the matrix equality Ax=B.
- _get_solution: To get the solution for the interpolation

In almost all cases, the user should not use this module.
The functionality of interpolation is nicely wrapped in another module, SplineInterpolation.

Created on Dec 28 17:28:17 2023
"""

__all__ = ['SPLINE', 'LinearSpline', 'QuadraticSpline', 'NaturalCubicSpline']

import numpy as np
from custom_inherit import doc_inherit

from .. import DOC_STYLE, FList, FListOrLList, IFloat, IFloatOrFList, LList, LLList, OptFunc, OptList, TOLERANCE
from ..__backend.errors_ import AtLeastOneParameterRequired, WrongBoundaryEquation


class SPLINE:
    """
    Base class for implementation of SPLINE interpolation method.

    Methods
    -------
    show_splines:
        Shows the spline equations.
    interpolate:
        Interpolates the ``value_to_approximate``.

    Notes
    -----
    Either ``function`` or ``function_values`` must be specified.
    """

    def __init__(self, given_values: FList, value_to_approximate: IFloat, function: OptFunc = None,
                 function_values: OptList = None):
        """

        Parameters
        ----------
        given_values:
            A list of x values.
        value_to_approximate:
            The value to approximate using interpolation method.
        function:
            A function to apply on the given values to get the resultant values.
        function_values:
            A list of values corresponding to the given values.

        Raises
        ------
        AtLeastOneParameterRequired:
            If both ``function`` and ``function_value`` parameters are not provided.

        ValueError:
            If ``value_to_approximate`` is not in the list of ``given_values``.
        """
        if function is None and function_values is None:
            raise AtLeastOneParameterRequired("One of ``function`` or ``function_values`` parameter is required.")

        if value_to_approximate < any(x for x in given_values) < value_to_approximate:
            raise ValueError("The given value to approximate is out of bounds from the given values.")

        self.given_values = given_values
        self.value_to_approximate = value_to_approximate
        self.function_values = function_values if function_values else [function(value) for value in given_values]

    def _solve(self):
        pass

    def show_splines(self):
        """Prints the spline equations."""

    def interpolate(self):
        """Interpolates the spline equations."""


@doc_inherit(SPLINE, style=DOC_STYLE)
class LinearSpline(SPLINE):
    """Implements the linear spline interpolation."""

    @doc_inherit(SPLINE.__init__, style=DOC_STYLE)
    def __init__(self, given_values: FList, value_to_approximate: IFloat, function: OptFunc = None,
                 function_values: OptList = None):
        """Initializes the linear spline interpolation class."""

        super().__init__(given_values, value_to_approximate, function, function_values)

    def _solve(self, give_values: bool = False) -> IFloatOrFList:
        """
        Solver method for linear interpolation.

        Parameters
        ----------
        give_values:
            Whether to return only the spline value or with the x/y values.
        """

        given, func_vals = self.given_values, self.function_values
        x_to_approx, n_points = self.value_to_approximate, len(self.given_values) - 1

        x_values = [(given[i], given[i + 1]) for i in range(n_points)]
        y_values = [(func_vals[i], func_vals[i + 1]) for i in range(n_points)]

        spline = [((x_to_approx - x2) * y1 - (x_to_approx - x1) * y2) / (x1 - x2)
                  for (x1, x2), (y1, y2) in zip(x_values, y_values)]

        return [x_values, y_values, spline] if give_values else spline

    def show_splines(self) -> None:
        x_values, y_values, spline = self._solve(give_values=True)

        print('The splines are approximated to 4 decimal places for display purposes only.')
        for i in range(len(spline)):
            den1 = x_values[i][0] - x_values[i][1]
            den2 = -den1

            temp1_ = (y_values[i][0] / den1) + (y_values[i][1] / den2)

            temp2_ = (y_values[i][0] * x_values[i][1]) / den1
            temp2_ += (y_values[i][1] * x_values[i][0]) / den2
            temp2_ *= -1

            print(f'Sp{i + 1}: {temp1_:+08.4f}x {temp2_:+08.4f}'
                  f'\t; x ∈ {"[" if i == 0 else "("}{self.given_values[i]:+}, {self.given_values[i + 1]:+}]')

    def interpolate(self) -> IFloat:
        """
        Performs the linear spline interpolations and return the interpolated value.

        Returns
        -------
            The interpolated value.
        """

        return _get_solution(self.given_values, self.value_to_approximate, self._solve())


@doc_inherit(SPLINE, style=DOC_STYLE)
class QuadraticSpline(SPLINE):
    """Implements the quadratic spline interpolation."""

    @doc_inherit(SPLINE.__init__, style=DOC_STYLE)
    def __init__(self, given_values, value_to_approximate, function=None, function_values=None, last_equation='first'):
        """Initializes the quadratic spline interpolation class."""

        super().__init__(given_values, value_to_approximate, function, function_values)
        self.last_equation = last_equation

    def _solve(self):
        given, fn_values = self.given_values, self.function_values
        n_given = len(given)
        len_mat = 3 * (n_given - 1)

        matrix_a = [[0] * len_mat for _ in range(len_mat)]
        matrix_b = [0] * len_mat

        i_index = _fill_initial_splines(given, fn_values, matrix_a, matrix_b, 0, 2)
        combined = _generate_middle_point_equations(given, 2)
        _fill_middle_splines(combined, i_index, matrix_a, n_given, 2)

        if self.last_equation == 'first':
            matrix_a[-1][0] = 1
        elif self.last_equation == 'last':
            matrix_a[-1][-3] = 1
        else:
            raise WrongBoundaryEquation('Possible values are \'first\' or \'last\'.')

        solution = _matrix_solver(matrix_a, matrix_b, n_given - 1, 2)

        return solution

    def show_splines(self):
        given, solution = self.given_values, self._solve()

        print('The splines are approximated to 4 decimal places for display purposes only.')
        for index, value in enumerate(solution):
            print(f'Sp{index + 1}: {value[0]:+08.4f}x^2 {value[1]:+08.4f}x {value[2]:+08.4f}'
                  f'\t; x ∈ {"[" if index == 0 else "("}{given[index]:+}, {given[index + 1]:+}]')

    def interpolate(self, n_derivative: IFloat = 0) -> IFloat:
        """
        Performs the quadratic spline interpolation and return the interpolated value.

        Parameters
        ----------
        n_derivative:
            The value of spline derivative to return. Default is 0, the original quadratic spline.

        Returns
        -------
            The interpolated value.
        """

        approx = self.value_to_approximate
        solution = _get_solution(self.given_values, approx, self._solve())

        if n_derivative == 0:
            return solution[0] * approx**2 + solution[1] * approx + solution[2]
        elif n_derivative == 1:
            return 2 * solution[0] * approx + solution[1]
        elif n_derivative == 2:
            return 2 * solution[0]
        else:
            return 0


@doc_inherit(SPLINE, style=DOC_STYLE)
class NaturalCubicSpline(SPLINE):
    """Implements the natural cubic spline interpolation."""

    @doc_inherit(SPLINE.__init__, style=DOC_STYLE)
    def __init__(self, given_values, value_to_approximate, function=None, function_values=None, last_equation='first'):
        """Initializes the natural cubic spline interpolation class."""

        super().__init__(given_values, value_to_approximate, function, function_values)
        self.last_equation = last_equation

    def _solve(self):
        given, fn_values = self.given_values, self.function_values
        n_given = len(given)
        len_mat = 4 * (n_given - 1)

        matrix_a = [[0] * len_mat for _ in range(len_mat)]
        matrix_b = [0] * len_mat

        i_index1 = _fill_initial_splines(given, fn_values, matrix_a, matrix_b, 0, 3)
        combined = _generate_middle_point_equations(given, 3)
        combined_der = _generate_middle_point_equations(given, -1)

        i_index1 = _fill_middle_splines(combined, i_index1, matrix_a, n_given, 3)
        _fill_middle_splines(combined_der, i_index1, matrix_a, n_given, 3)

        matrix_a[-2][0] = 6 * matrix_a[0][2]
        matrix_a[-2][1] = 2 * matrix_a[0][3]
        matrix_a[-1][-4] = 6 * matrix_a[2 * (n_given - 1) - 1][-2]
        matrix_a[-1][-3] = 2 * matrix_a[2 * (n_given - 1) - 1][-1]

        solution = _matrix_solver(matrix_a, matrix_b, n_given - 1, 3)

        return solution

    def show_splines(self):
        given, solution = self.given_values, self._solve()

        print('The splines are approximated to 4 decimal places for display purposes only.')
        for index, value in enumerate(solution):
            print(f'Sp{index + 1}: {value[0]:+08.4f}x^3 {value[1]:+08.4f}x^2 {value[2]:+08.4f}x '
                  f'{value[3]:+08.4f}\t; x ∈ {"[" if index == 0 else "("}{given[index]:+}, {given[index + 1]:+}]')

    def interpolate(self, n_derivative: IFloat = 0) -> IFloat:
        """
        Performs the natural cubic spline interpolations and return the interpolated value.

        Parameters
        ----------
        n_derivative:
            The value of spline derivative to return. Default is 0, the original cubic spline.

        Returns
        -------
            The interpolated value.
        """

        approx = self.value_to_approximate
        solution = _get_solution(self.given_values, approx, self._solve())

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


def _fill_initial_splines(given_values: FList, function_values: FList, matrix_a: LList, matrix_b: FList,
                          i_index: IFloat, deg_spline: IFloat) -> IFloat:
    """
    Fills in the initial spline equations in the matrix.

    Parameters
    ----------
    given_values:
        A list of given x values.
    function_values:
        A list of corresponding function values against the x values.
    matrix_a:
        The matrix A, in Ax=B.
    matrix_b:
        The matrix B, in Ax=B.
    i_index:
        The current iterative index. The index is updated and returned to be used further.
    deg_spline:
        The degree of the spline. 2 for quadratic and 3 for cubic spline.

    Returns
    -------
        Updated value of i_index1.
    """

    condition, deg_spline = 2 * (len(given_values) - 1), deg_spline + 1
    j_index, k_index = 0, 0

    while i_index < condition:
        ind = deg_spline * k_index
        val1, val2 = given_values[j_index], given_values[j_index + 1]

        if deg_spline == 3:
            matrix_a[i_index][ind:ind + deg_spline] = [val1**2, val1, 1]
            matrix_a[i_index + 1][ind:ind + deg_spline] = [val2**2, val2, 1]
        else:
            matrix_a[i_index][ind:ind + deg_spline] = [val1**3, val1**2, val1, 1]
            matrix_a[i_index + 1][ind:ind + deg_spline] = [val2**3, val2**2, val2, 1]

        matrix_b[i_index] = function_values[j_index]
        matrix_b[i_index + 1] = function_values[j_index + 1]

        i_index += 2
        j_index += 1
        k_index += 1

    return i_index


def _generate_middle_point_equations(given_values: FList, deg_spline: IFloat) -> LLList:
    """
    Creates the middle spline equation for quadratic and cubic splines.

    Parameters
    ----------
    given_values:
        The provided list of x values.
    deg_spline:
        The degree of the spline. 2 for quadratic and 3 for cubic spline.

    Returns
    -------
        A list containing all the second (or third, for cubic spline) derivative values of quadratic and cubic splines.
    """

    given, deg_spline = given_values[1:-1], deg_spline + 1

    if deg_spline == 3:
        a_coefficients = [[2 * j * i for i in given] for j in [1, -1]]
        b_coefficients = [[j for _ in given] for j in [1, -1]]
        c_coefficients = [[0] * len(given) for _ in range(2)]

        res = [a_coefficients, b_coefficients, c_coefficients]
    elif deg_spline == 4:
        a_coefficients = [[3 * j * i**2 for i in given] for j in [1, -1]]
        b_coefficients = [[2 * j * _ for _ in given] for j in [1, -1]]
        c_coefficients = [[j for _ in given] for j in [1, -1]]
        d_coefficients = [[0] * len(given) for _ in range(2)]

        res = [a_coefficients, b_coefficients, c_coefficients, d_coefficients]
    else:
        a_coefficients = [[6 * j * i for i in given] for j in [1, -1]]
        b_coefficients = [[2 * j for _ in given] for j in [1, -1]]
        c_coefficients = [[0] * len(given) for _ in range(2)]
        d_coefficients = c_coefficients

        res = [a_coefficients, b_coefficients, c_coefficients, d_coefficients]

    return [[[*p] for p in zip(*_)] for _ in zip(*res)]


def _fill_middle_splines(combined_values: LLList, i_index1: IFloat, matrix_a: LList, n_given: IFloat,
                         deg_spline: IFloat) -> IFloat:
    """
    Fills in the middle splines equations in the matrix.

    Parameters
    ----------
    combined_values:
        A list containing sub-lists with spline coefficients.
    i_index1:
        Starting index to fill in the matrix.
    matrix_a:
        Matrix to be filled with spline coefficients.
    n_given:
        Number of given values.
    deg_spline:
        The degree of the spline. 2 for quadratic and 3 for cubic spline.

    Returns
    -------
        Updated value of i_index1.
    """

    deg_spline, j_index, k_index = deg_spline + 1, 0, 0
    i_temp = i_index1

    while i_index1 < i_temp + n_given - 2:
        ind = deg_spline * k_index
        matrix_a[i_index1][ind:ind + deg_spline] = combined_values[0][j_index]
        matrix_a[i_index1][ind + deg_spline:ind + 2 * deg_spline] = combined_values[1][j_index]

        i_index1 += 1
        j_index += 1
        k_index += 1

    return i_index1


def _matrix_solver(matrix_a: LList, matrix_b: FList, n_given: IFloat, deg_spline: IFloat):
    """
    Solves the matrix for the interpolation method.

    Parameters
    ----------
    matrix_a:
        The matrix A, in Ax=B.
    matrix_b:
        The matrix B, in Ax=B.
    n_given:
        Number of given data points.
    deg_spline:
        The degree of the spline. 2 for quadratic and 3 for cubic spline.

    Returns
    -------
        A list of values containing the result of matrix equation Ax=B.
    """

    # solve using numpy solver
    solution = np.linalg.solve(np.array(matrix_a), np.array(matrix_b))
    # reshape according to the data
    solution = np.reshape(solution, (n_given, deg_spline + 1))
    # convert all very small values to 0
    solution[np.abs(solution) < TOLERANCE] = 0

    return solution.tolist()


def _get_solution(given_values: FList, value_to_approximate: IFloat, solution: FListOrLList) -> IFloatOrFList:
    """
    Get the interpolated/approximated value for the specified ``value_to_approximate`` using the given solution.

    Parameters
    ----------
    given_values:
        A list of given x values.
    value_to_approximate:
        The x value for which interpolation/approximation is required.
    solution:
        The solution obtained from the interpolation method.

    Returns
    -------
        The interpolated/approximated value for the specified `value_to_approximate`.
    """

    given_values.insert(0, value_to_approximate)
    given_values.sort()
    index_to_remove = given_values.index(value_to_approximate) - 1
    del given_values[index_to_remove + 1]

    return solution[index_to_remove]
