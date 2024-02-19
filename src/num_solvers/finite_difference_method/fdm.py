"""Created on Jan 11 13:59:35 2024"""

from umatrix.matrix import null_matrix

from src.num_solvers import FList, IFloat, IFloatOrFList


def tri_diagonal_matrix(elements, n_rows, n_cols, diff):
    zero = null_matrix(n_rows, n_cols)

    for row in range(zero.n_rows):
        for col in range(zero.n_cols):
            if row == col + diff:
                zero[row][col] = elements[0]
            if row == col:
                zero[row][col] = elements[1]
            if row == col - diff:
                zero[row][col] = elements[2]

    return zero


def tri_diagonal_matrix2(elements, n_rows, n_cols):
    zero = null_matrix(n_rows, n_cols)

    for row in range(zero.n_rows - 1):
        for i, v in enumerate(elements):
            zero[row][i + row - 1] = v

    for i, v in enumerate(elements[:-1]):
        zero[-1][-i - 1] = v

    zero[0][-1] = 0

    return zero


def d2_matrix(time_step: int, space_step: int, d_type: str = 'x', delta_factor: IFloat = 1):
    values = [i / delta_factor**2 for i in [1, -2, 1]]
    return tri_diagonal_matrix(values, time_step, space_step, 1 if d_type == 'x' else 4)


def d1_central(time_step: int, space_step: int, d_type: str = 'x', delta_factor: IFloat = 1):
    values = [i / (2 * delta_factor) for i in [-1, 0, 1]]
    return tri_diagonal_matrix(values, time_step, space_step, 1 if d_type == 'x' else 4)


def bi_diagonal_space_matrix(n_rows, n_cols):
    zero = null_matrix(n_rows, n_cols)

    for row in range(zero.n_rows):
        for col in range(zero.n_cols):
            if row == col - 1:
                zero[row][col] = 1
            if row == col:
                zero[row][col] = -1

    return zero


def d1_fwd(time_step: int, space_step: int, d_type: str = 'x', delta_factor: IFloat = 1):
    values = [i / delta_factor for i in [-1, 0, 1]]
    return tri_diagonal_matrix(values, time_step, space_step, 1 if d_type == 'x' else 4)


def d1_bkw(time_step: int, space_step: int, d_type: str = 'x', delta_factor: IFloat = 1):
    values = [i / delta_factor for i in [1, 0, -1]]
    return tri_diagonal_matrix(values, time_step, space_step, 1 if d_type == 'x' else 4)


def backwards_euler(time_step: int, space_step: int, d_type: str = 'x', delta_factor: IFloat = 1):
    values = [i / delta_factor for i in [0, -1, 1]]
    return tri_diagonal_matrix(values, time_step, space_step, 1 if d_type == 'x' else 4)


def single_term(time_step: int, space_step: int, delta_factor: IFloat = 1):
    return tri_diagonal_matrix([0, 1, 0], time_step, space_step, 1) * delta_factor


def boundary_matrix(time_step: int, boundary_condition: FList):
    zeros = null_matrix(time_step, 1)

    zeros[0] = [boundary_condition[0]]
    zeros[-1] = [boundary_condition[1]]

    return zeros


def apply_ic(time_step: int, initial_condition: IFloatOrFList):
    zeros = null_matrix(time_step, 1).t

    if isinstance(initial_condition, float):
        for i in range(zeros.n_cols):
            zeros[i] = initial_condition
    else:
        zeros.elements = initial_condition

    return zeros.t


class SingleVariableODESolver:

    def __init__(self, fdm_functions, delta, steps, boundary_conditions, has_single_term=True):
        self.fdm = fdm_functions
        self.delta = delta
        self.steps = steps
        self.bc = boundary_conditions
        self.hST = has_single_term

    def lhs(self):
        if self.hST:
            temp_ = ([i(self.steps, self.steps, self.delta) for i in self.fdm[:-1]] +
                     [self.fdm[-1](self.steps, self.steps, 1)])
        else:
            temp_ = [i(self.steps, self.steps, self.delta) for i in self.fdm[:-1]]

        return temp_

    def rhs(self):
        return boundary_matrix(self.steps, 1)
