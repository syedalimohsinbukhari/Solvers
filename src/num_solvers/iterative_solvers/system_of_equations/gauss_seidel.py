"""Gauss-Seidel

This module provides functionality to solve diagonally dominant system of equations using Gauss-Seidel method.

Created on Oct 14 05:46:04 2023
"""

__all__ = ['gauss_seidel']

from custom_inherit import doc_inherit

from ... import DOC_STYLE, FList, FTuple, IFloat, LLList, LList, TOLERANCE
from ...__backend.iterative_solvers_ import GaussSeidel, SysEqnSolver


@doc_inherit(SysEqnSolver.__init__, style=DOC_STYLE)
def gauss_seidel(system_of_equations: LList,
                 solution: FList,
                 n_iter: int = 500,
                 tolerance: IFloat = TOLERANCE,
                 initial_guess: FTuple = (0, 0, 0),
                 n_decimal: int = 8,
                 get_solution_set: bool = False) -> LLList:
    """
    Solve a system of linear equations using the Gauss-Seidel method.

    Parameters
    ----------

    get_solution_set:
        If True, returns both the solution and the solution set for each iteration, default is False.

    Returns
    -------
        If ``get_solution_set`` is ``False``, returns the solution as a list.
        Otherwise, returns a tuple containing the solution and a list of solution sets for each iteration.

    Examples
    --------
    >>> equations = [[2, 1, 1], [3, 10, 2], [2, 1, 4]]
    >>> solutions = [5, 10, 9]
    >>> GS_Solver = gauss_seidel(equations, solutions)
    """

    gs = GaussSeidel(system_of_equations, solution, n_iter, tolerance, initial_guess, n_decimal)

    return gs.solve() if not get_solution_set else (gs.solve(), gs.solution_set)
