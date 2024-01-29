"""Iterative Solvers

This module provides functionality for iterative matrix solvers, such as,

- GaussJacobi
- GaussSeidel

using the base class,

- SysEqnSolver

Currently, this module supports only 3x3, diagonally dominant matrices.

Created on Oct 14 06:07:34 2023
"""

__all__ = ['SysEqnSolver', 'GaussJacobi', 'GaussSeidel']

from .core_helpers_ import round_list_
from .. import FList, FTuple, IFloat, LLList, LList, N_DECIMAL, TOLERANCE


# TODO: re-implement SysEqnSolver using Matrix class internally.


class SysEqnSolver:
    """Base class for system of equation solvers."""

    def __init__(self, system_of_equations: LList, solution: FList, n_iter: int = 500, tolerance: IFloat = TOLERANCE,
                 initial_guess: FTuple = (0, 0, 0), n_decimal: int = N_DECIMAL):
        """
        Initializer for SysEqnSolver class.

        Parameters
        ----------
        system_of_equations:
            A list of system of equations.
        solution:
            The solution to the provided system of equations.
        n_iter:
            Number of iterations to perform. Defaults to 500.
        tolerance:
            The tolerance level for solution. Defaults to 1e-8.
        initial_guess:
            The initial guess for the solutions to the system of equations. Defaults to (0, 0, 0).
        n_decimal:
            Number of digits to round off to. Default is 8.
        """

        self.n_iter = n_iter
        self.solution = solution
        self.n_decimal = n_decimal
        self.tolerance = tolerance
        self.initial_guess = initial_guess
        self.system_of_equations = system_of_equations
        self._arr1, self._arr2, self._arr3 = [self.initial_guess[0]], [self.initial_guess[1]], [self.initial_guess[2]]

        self.__sanity_check()

    def __sanity_check(self):
        def _auxiliary(val1, val2, val3):
            return abs(val1) >= abs(val2) + abs(val3)

        sys_eq = self.system_of_equations
        cond1 = _auxiliary(sys_eq[0][0], sys_eq[0][1], sys_eq[0][2])
        cond2 = _auxiliary(sys_eq[1][1], sys_eq[1][0], sys_eq[1][2])
        cond3 = _auxiliary(sys_eq[2][2], sys_eq[2][0], sys_eq[2][1])

        if cond1 and cond2 and cond3:
            pass
        else:
            print('The system of equations is not diagonally dominant. The solutions might not be correct.')

    @property
    def solution_set(self) -> LList:
        """
        Gives the solution set, e.g., all the iterations for the system of equations.

        Returns
        -------
            List of all iterations for the solution.
        """

        if not self._arr1 or len(self._arr1) == 1:
            self.solve()

        return [self._arr1, self._arr2, self._arr3]

    def _evaluate(self):
        class_name, tolerance = self.__class__.__name__, self.tolerance
        equations, solution, initial_guess = self.system_of_equations, self.solution, self.initial_guess

        iter1_ = solution[0]
        iter1_ -= equations[0][1] * initial_guess[1]
        iter1_ -= equations[0][2] * initial_guess[2]
        iter1_ = iter1_ * equations[0][0]**-1

        iter2_ = solution[1]

        if class_name == 'GaussJacobi':
            iter2_ -= equations[1][0] * initial_guess[0]
        else:
            iter2_ -= equations[1][0] * iter1_

        iter2_ -= equations[1][2] * initial_guess[2]
        iter2_ = iter2_ * equations[1][1]**-1

        iter3_ = solution[2]

        if class_name == 'GaussJacobi':
            iter3_ -= equations[2][0] * initial_guess[0]
            iter3_ -= equations[2][1] * initial_guess[1]
        else:
            iter3_ -= equations[2][0] * iter1_
            iter3_ -= equations[2][1] * iter2_

        iter3_ = iter3_ * equations[2][2]**-1

        self._arr1.append(iter1_)
        self._arr2.append(iter2_)
        self._arr3.append(iter3_)

        return [iter1_, iter2_, iter3_]

    def __break_iteration(self):
        cond1 = abs(self._arr1[-1] - self._arr1[-2]) < self.tolerance
        cond2 = abs(self._arr2[-1] - self._arr2[-2]) < self.tolerance
        cond3 = abs(self._arr3[-1] - self._arr3[-2]) < self.tolerance

        return cond1 and cond2 and cond3

    def solve(self) -> LLList:
        """
        Solves the system of equations.

        Returns
        -------
            Solutions for the provided system of equation.
        """

        n_iter, n_decimal = self.n_iter, self.n_decimal

        for iter_ in range(n_iter):
            self.initial_guess = self._evaluate()
            if self.__break_iteration():
                break

        return round_list_(self.initial_guess, n_decimal)


class GaussJacobi(SysEqnSolver):
    """This class implements the Gauss-Jacobi method of solving system of equations."""
    pass


class GaussSeidel(SysEqnSolver):
    """This class implements the Gauss-Seidel method of solving system of equations."""
    pass
