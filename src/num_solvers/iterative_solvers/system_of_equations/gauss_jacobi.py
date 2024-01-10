"""Gauss-Jacobi

This module provides functionality to solve diagonally dominant system of equations using Gauss-Jacobi method.

Created on Oct 13 11:37:27 2023
"""

__all__ = ['gauss_jacobi']

from typing import List, Tuple, Union

import numpy as np
from custom_inherit import doc_inherit

from ... import DOC_STYLE, FList, FTuple, IFloat, LLList, LList, TOLERANCE
from ...__backend.iterative_solvers_ import GaussJacobi, SysEqnSolver

GJ_OUTPUT = Union[np.ndarray, Tuple[np.ndarray, List[np.ndarray]]]
INITIAL_GUESS = Tuple[float, float, float]


@doc_inherit(SysEqnSolver.__init__, style=DOC_STYLE)
def gauss_jacobi(system_of_equations: LList,
                 solution: FList,
                 n_iter: int = 500,
                 tolerance: IFloat = TOLERANCE,
                 initial_guess: FTuple = (0, 0, 0),
                 n_decimal: int = 8,
                 get_solution_set: bool = False) -> LLList:
    """
    Solve a system of linear equations using the Gauss-Jacobi method.

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
    >>> GJ_Solver = gauss_jacobi(equations, solutions)
    """

    gj = GaussJacobi(system_of_equations, solution, n_iter, tolerance, initial_guess, n_decimal)

    return gj.solve() if not get_solution_set else (gj.solve(), gj.solution_set)
