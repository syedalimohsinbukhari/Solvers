"""Created on Jan 22 23:00:22 2024"""

from src.num_solvers.matrix_decomposition.matrix import Matrix


def steepest_descent(matrix: Matrix, solution: Matrix, initial_guess: Matrix, n_iterations: int = 50):
    r = [0] * n_iterations
    new_guess = [0] * n_iterations

    new_guess[0] = guess

    r[0] = solution - matrix * guess[0]
    print(r[0])


matrix = Matrix([[-10, -4, 14], [-2, -4, 6]])
solution = Matrix([[2], [-8]])
guess = Matrix([0, 0])

steepest_descent(matrix, solution, guess)
