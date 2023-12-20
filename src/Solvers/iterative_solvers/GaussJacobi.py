"""Created on Oct 13 11:37:27 2023"""

from typing import List

from .backend.ITERATIVE_SOLVER_ import GaussJacobi


def gauss_jacobi(system_of_equations: List[List], solution: List, n_iter: int = 500, tol: float = 1e-5,
                 initial_guess: tuple = (0, 0, 0), get_solution_set=False):
    """

    Examples
    --------
    equations = [[2, 1, 1], [3, 10, 2], [2, 1, 4]]
    solutions = [5, 10, 9]

    GJ_Solver = gauss_jacobi(equations, solutions)
    """
    gj = GaussJacobi(system_of_equations, solution, n_iter, tol, initial_guess)

    return gj.solve() if not get_solution_set else (gj.solve(), gj.solution_set)
