"""Created on Oct 14 05:46:04 2023"""

from typing import List

from .backend.ITERATIVE_SOLVER_ import GaussSeidel


def gauss_seidel(system_of_equations: List[List], solution: List, n_iter: int = 500, tol: float = 1e-5,
                 initial_guess: tuple = (0, 0, 0), get_solution_set=False):
    """

    Examples
    --------
    equations = [[2, 1, 1], [3, 10, 2], [2, 1, 4]]
    solutions = [5, 10, 9]

    GJ_Solver = gauss_seidel(equations, solutions)
    """
    gs = GaussSeidel(system_of_equations, solution, n_iter, tol, initial_guess)

    return gs.solve() if not get_solution_set else (gs.solve(), gs.solution_set)
