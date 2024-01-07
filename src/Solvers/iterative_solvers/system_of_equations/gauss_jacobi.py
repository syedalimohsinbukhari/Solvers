"""Created on Oct 13 11:37:27 2023"""

from typing import List, Tuple, Union

import numpy as np

from ...__backend.iterative_solvers_ import GaussJacobi

GJ_OUTPUT = Union[np.ndarray, Tuple[np.ndarray, List[np.ndarray]]]
INITIAL_GUESS = Tuple[float, float, float]


def gauss_jacobi(system_of_equations: List[List[float]],
                 solution: List[float],
                 n_iter: int = 500,
                 tol: float = 1e-5,
                 initial_guess: INITIAL_GUESS = (0, 0, 0),
                 get_solution_set: bool = False) -> GJ_OUTPUT:
    """
    Solve a system of linear equations using the Gauss-Jacobi method.

    Parameters
    ----------
    system_of_equations : List[List[float]]
        Coefficient matrix of the system of linear equations.
    solution : List[float]
        Right-hand side of the system of linear equations.
    n_iter : int, optional
        Maximum number of iterations, default is 500.
    tol : float, optional
        Tolerance for convergence, default is 1e-5.
    initial_guess : INITIAL_GUESS, optional
        Initial guess for the solution, default is (0, 0, 0).
    get_solution_set : bool, optional
        If True, returns both the solution and the solution set for each iteration, default is False.

    Returns
    -------
    GJ_Output
        If ``get_solution_set`` is ``False``, returns the solution as a NumPy array.
        If ``get_solution_set`` is ``True``, returns a tuple containing the solution and a list of solution sets for
        each iteration.

    Examples
    --------
    >>> equations = [[2, 1, 1], [3, 10, 2], [2, 1, 4]]
    >>> solutions = [5, 10, 9]
    >>> GJ_Solver = gauss_jacobi(equations, solutions)
    """
    gj = GaussJacobi(system_of_equations, solution, n_iter, tol, initial_guess)

    return gj.solve() if not get_solution_set else (gj.solve(), gj.solution_set)
