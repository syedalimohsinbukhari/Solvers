"""Created on Feb 19 23:52:28 2024"""

__all__ = ['OneDimensionalFDM', 'OneDimensionalPDESolver', 'bi_diagonal_matrix', 'enforce_boundary_condition',
           'initial_condition_matrix', 'tri_diagonal_matrix', 'mesh2d']

from math import floor

from umatrix.matrix import identity_matrix, Matrix, null_matrix

from .core_helpers_ import round_value_
from .matrix_ import matrix_copy
from .. import FList, Func, IFloat, IFloatOrFList, LMat, N_DECIMAL, OptIFloat, OptList, TOLERANCE


class OneDimensionalFDM:

    def __init__(self, x_range, delta_x, delta_t, n_steps: OptIFloat = None, wrap_boundaries: bool = False,
                 elements: OptList = None, tolerance: IFloat = TOLERANCE):
        self.x_range = x_range
        self.dx = delta_x
        self.dt = delta_t
        self.wrap_boundaries = wrap_boundaries
        self.elements = elements

        # take the difference of provided x values
        x_diff = x_range[1] - x_range[0]

        # Calculate the number of steps based on the specified step size or the provided x range
        self.n_steps = n_steps if n_steps is not None else floor(x_diff / delta_x)

        # Calculate the actual step size based on the calculated number of steps
        dx2 = x_diff / self.n_steps

        # Adjust the step size if it differs from the specified step size
        if abs(delta_x - dx2) > tolerance:
            print('The recalculated `delta_x` from `n_steps` does not match the provided value.\n'
                  f'The `delta_x` parameter has been modified from {delta_x} to '
                  f'{round_value_(dx2, N_DECIMAL)}')
            self.dx = dx2

    def __factors(self, diff_type):
        dx, dt = self.dx, self.dt
        constants = {'fwd': dt / dx,
                     'bkw': dt / dx,
                     'cnt': dt / (2 * dx),
                     'cnt2': dt / dx**2}

        return constants[diff_type]

    def d1_forward(self):
        return self.__factors('fwd') * bi_diagonal_matrix(self.n_steps,
                                                          self.wrap_boundaries,
                                                          'fwd',
                                                          self.elements)

    def d1_backward(self):
        return self.__factors('bkw') * bi_diagonal_matrix(self.n_steps,
                                                          self.wrap_boundaries,
                                                          'bkw',
                                                          self.elements)

    def d1_central(self):
        return self.__factors('cnt') * bi_diagonal_matrix(self.n_steps,
                                                          self.wrap_boundaries,
                                                          'cnt',
                                                          self.elements)

    def d2_central(self):
        return self.__factors('cnt2') * tri_diagonal_matrix(self.n_steps,
                                                            self.wrap_boundaries,
                                                            self.elements)


class OneDimensionalPDESolver:

    def __init__(self, initial_condition, boundary_conditions, fdm_matrices, has_single_term=True,
                 ic_values: OptList = None):
        self.ic = initial_condition
        self.bc = boundary_conditions
        self.fdm = fdm_matrices
        self.hST = has_single_term
        self.ic_values = ic_values

        self.n_steps = self.fdm[0].n_cols

    def lhs(self):
        identity_ = identity_matrix(self.n_steps)
        for matrix_ in self.fdm[1:]:
            identity_ += matrix_

        return identity_

    def rhs(self):
        return initial_condition_matrix(self.n_steps, self.ic, self.ic_values)

    def __solver(self):
        return self.lhs().inverse() * self.rhs()

    def solve(self, time_steps: OptIFloat = 10):
        lhs, temp_ = self.lhs().inverse(), self.__solver()

        solution = [temp_]

        for i in range(1, time_steps):
            solution.append(lhs * solution[i - 1])
            enforce_boundary_condition(solution[i], self.bc, i > 0)

        return solution


def initial_condition_matrix(n_steps: int, initial_condition: IFloatOrFList or Func, values: OptList = None):
    ic_matrix = null_matrix(n_steps, 1).t

    if isinstance(initial_condition, (int, float)):
        ic_matrix.elements = [initial_condition] * n_steps

    elif isinstance(initial_condition, list):
        ic_matrix.elements = initial_condition

    elif isinstance(initial_condition, Func):
        if len(values) != n_steps:
            raise ValueError('The length of vector provided does not match with the number of steps provided')
        ic_matrix.elements = [initial_condition(i) for i in values]

    return Matrix(ic_matrix).t


def enforce_boundary_condition(matrix: Matrix, boundary_conditions: FList, overwrite: bool = False):
    matrix_ = matrix_copy(matrix, overwrite) if overwrite else matrix
    if overwrite:
        matrix_[0][0] = boundary_conditions[0]
        matrix_[-1][0] = boundary_conditions[-1]

    return matrix_


def bi_diagonal_matrix(n_steps, wrap_boundaries: bool = False, diff_type: str = 'fwd',
                       elements: OptList = None) -> Matrix:
    difference_indices = {'fwd': [-1, 0],
                          'bkw': [1, 0],
                          'cnt': [-1, 1]}

    difference_ = difference_indices[diff_type]

    elements = elements if elements else [1, -1]
    bi_diagonal_ = null_matrix(n_steps)

    for row in range(bi_diagonal_.n_rows):
        for col in range(bi_diagonal_.n_cols):
            if row == col + difference_[0]:
                bi_diagonal_[row][col] = elements[0]
            if row == col + difference_[1]:
                bi_diagonal_[row][col] = elements[1]

    if wrap_boundaries:
        if diff_type == 'bkw':
            bi_diagonal_[0][-1] = elements[0]
        elif diff_type == 'fwd':
            bi_diagonal_[-1][0] = elements[0]
        else:
            bi_diagonal_[0][-1] = elements[1]
            bi_diagonal_[-1][0] = elements[0]

    return bi_diagonal_


def tri_diagonal_matrix(n_steps, wrap_boundaries: bool = False, elements: OptList = None) -> Matrix:
    elements = elements if elements else [1, -2, 1]
    tri_diagonal_ = null_matrix(n_steps)

    for row in range(tri_diagonal_.n_rows):
        for col in range(tri_diagonal_.n_cols):
            if row == col + 1:
                tri_diagonal_[row][col] = elements[0]
            if row == col:
                tri_diagonal_[row][col] = elements[1]
            if row == col - 1:
                tri_diagonal_[row][col] = elements[2]

    if wrap_boundaries:
        tri_diagonal_[0][-1] = elements[-1]
        tri_diagonal_[-1][0] = elements[-1]

    return tri_diagonal_


def mesh2d(x: FList, y: FList) -> LMat:
    """Create a meshgrid from two 1D arrays."""
    if isinstance(x, Matrix):
        x = x.elements

    if isinstance(y, Matrix):
        y = y.elements

    x_grid = [[xi for _ in range(len(y))] for xi in x]
    y_grid = [[yj for yj in y] for _ in range(len(x))]

    return [Matrix(x_grid).t, Matrix(y_grid).t]
