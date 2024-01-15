"""Created on Oct 04 21:36:44 2023"""

from fractions import Fraction
from math import sqrt

from .. import LList
from ..__backend import errors_ as err_


# TODO: Check the setting of values inside column matrices, they're acting up

class Matrix:

    def __init__(self, elements, n_decimal=-1):
        self.elements = elements
        self.n_decimal = n_decimal

    def __repr__(self):
        elements = self.elements

        out = ''
        if self.n_rows == 1:
            max_width = max(len(str(row)) for row in elements)
            out += '[' + ' '.join([str(element).rjust(max_width) for element in elements]) + ']'
        else:
            out = '['
            max_width = max(len(str(element)) for row in elements for element in row)
            for index, row in enumerate(elements):
                br = '[' if index == 0 else ' ['
                out += br + ' '.join([str(element).rjust(max_width) for element in row]) + ']'
                if index < self.n_rows - 1:
                    out += '\n'
            out += ']'

        return out

    def __str__(self):
        return repr(self)

    def __len__(self):
        return len(self.elements)

    def __eq__(self, other):
        return self.elements == other.elements

    def __ne__(self, other):
        return not self.__eq__(other)

    def __add__(self, other):
        if isinstance(other, Matrix):
            cond = self.n_rows == other.n_rows
            cond = cond and self.n_cols == other.n_cols

            if not cond:
                raise err_.MatrixDimensionsMismatch()

            if self._multi_rows() and other._multi_rows():
                row = [[self_element + other_element
                        for self_element, other_element in zip(self_row, other_row)]
                       for self_row, other_row in zip(self.elements, other.elements)]
            else:
                row = [self_element + other_element
                       for self_element, other_element in zip(self.elements, other.elements)]

            return self.give_output(row)

        elif isinstance(other, (int, float)):
            if self._multi_rows():
                result_elements = [[element + other for element in row] for row in self.elements]
            else:
                result_elements = [element + other for element in self.elements]

            return self.give_output(result_elements)

        else:
            raise ValueError("Unsupported operand type for addition.")

    def __sub__(self, other):
        return self.__add__(other=-1 * other)

    def __mul__(self, other):
        if isinstance(other, (int, float)):
            return self._scalar_vector_multiplication(other=other)

        elif isinstance(other, Matrix):
            if self.n_cols != other.n_rows:
                raise err_.MatrixDimensionsMismatch(f'Inner CxR={self.n_cols}x{other.n_rows}, not allowed.')

            if self._multi_rows():
                if other.n_cols == 1:
                    return self._multi_row_v_col_matrix(other)
                else:
                    return self._multi_matrix(other=other)
            else:
                if other._multi_rows():
                    return self._row_v_multi_row(other)
                else:
                    return self._row_v_col(other)

    __radd__ = __add__

    __rsub__ = __sub__

    __rmul__ = __mul__

    def __truediv__(self, other):
        if isinstance(other, Matrix):
            raise err_.DivisionByMatrix()

        if self._multi_rows():
            answer = [[element / other for element in row] for row in self.elements]
        else:
            answer = [element / other for element in self.elements]

        return self.give_output(answer)

    def __neg__(self):
        negated_elements = [[-element for element in row] for row in self.elements]
        return self.give_output(negated_elements)

    def __getitem__(self, index):
        if index > len(self.elements):
            raise err_.IndexOutOfBounds()

        if isinstance(index, slice):
            raise err_.SlicingNotAllowed()

        return self.give_output(self.elements[index]) if self._multi_rows() else self.elements[index]

    def __setitem__(self, index, value):
        if index >= len(self.elements):
            raise err_.IndexOutOfBounds()

        if isinstance(index, slice):
            raise err_.SlicingNotAllowed()

        self.elements[index] = value

    @property
    def n_rows(self):
        return len(self.elements) if self._multi_rows() else 1

    @property
    def n_cols(self):
        return len(self.elements[0]) if self._multi_rows() else len(self.elements)

    @property
    def dim(self):
        return f'RxC: {self.n_rows}x{self.n_cols}'

    @property
    def is_square(self):
        return self.n_rows == self.n_cols

    @property
    def is_singular(self):
        return True if determinant(self.elements) == 0 else False

    @property
    def trace(self):
        return sum([self.elements[i][i] for i in range(self.n_rows)])

    @property
    def in_fractions(self):
        fr_ = [[FRACTION(j).fraction for j in i] for i in self.elements]

        return self.give_output(fr_)

    @property
    def determinant(self):
        return determinant(self.elements)

    @property
    def inverse(self):
        inv_ = INVERSE(self.elements).inverse()
        return self.give_output(inv_)

    @property
    def t(self):
        return self.transpose()

    @property
    def adjoint_matrix(self):
        return (self.inverse * self.determinant).in_fractions

    @property
    def diagonal_of_matrix(self):
        n_rows, n_cols, elements = self.n_rows, self.n_cols, self.elements

        identity_ = identity_matrix(n_rows, n_cols)
        for i in range(n_rows):
            identity_[i][i] = elements[i][i]

        return identity_

    @staticmethod
    def give_output(output):
        return Matrix(output)

    def _multi_rows(self):
        return isinstance(self.elements[0], list)

    def _scalar_vector_multiplication(self, other):
        if self._multi_rows():
            result_elements = [[element * other for element in row] for row in self.elements]
        else:
            result_elements = [element * other for element in self.elements]

        return self.give_output(result_elements)

    def _row_v_col(self, other):
        row_v_col = [sum([k * l[0] for k, l in zip(self.elements, other.elements)])]

        return self.give_output(row_v_col)

    def _multi_matrix(self, other):
        result = []
        for row in self.elements:
            r_temp = []
            for c_row in other.t.elements:
                r_temp.append(sum([r_elem * c_elem for r_elem, c_elem in zip(row, c_row)]))
            result.append(r_temp)

        return self.give_output(result)

    def _row_v_multi_row(self, other):
        if other.n_cols == 1:
            mul_ = sum([sum([s_elem * o_elem[0]]) for o_elem, s_elem in zip(other.elements, self.elements)])
        else:
            temp_ = []
            for c_rw, r_rw in zip(self.elements, other.elements):
                temp_.append([c_rw * i for i in r_rw])

            mul_ = [sum([i[j] for i in temp_]) for j in range(len(temp_))]

        return self.give_output(mul_) if isinstance(mul_, list) and len(mul_) > 1 else mul_

    def _multi_row_v_col_matrix(self, other):
        multi_row_v_col_matrix = [[sum(self.elements[row][col] * other.elements[col][0]
                                       for col in range(self.n_cols))] for row in range(self.n_rows)]
        return self.give_output(multi_row_v_col_matrix)

    def _transpose(self):
        n_rows, n_cols, elements, give_output = self.n_rows, self.n_cols, self.elements, self.give_output

        if len(self) == 0:
            return give_output([])

        if self._multi_rows():
            answer = [[elements[row][col] for row in range(n_rows)] for col in range(n_cols)]

            if n_cols == 1:
                answer = answer[0]

            transposed_elements = answer
        else:
            if self.n_cols > 1:
                transposed_elements = [[element] for element in elements]
            else:
                transposed_elements = [element for element in elements]

        return give_output(transposed_elements)

    def transpose(self):
        return self._transpose()

    def hadamard_product(self, other):
        elements, o_elements = self.elements, other.elements

        result = [[i * j for i, j in zip(self_rows, other_rows)] for self_rows, other_rows in
                  zip(elements, o_elements)]

        return self.give_output(result)

    def elementwise_product(self, other):
        return self.hadamard_product(other)

    def is_multiplicative_inverse_of(self, other):
        return True if self * other == identity_matrix(self.n_rows, self.n_cols) else False

    def is_symmetric(self):
        return True if self == self.t else False

    def is_orthogonal(self):
        return True if self * self.t == identity_matrix(self.n_rows, self.n_cols) else False

    def is_orthogonal_to(self, other):
        if self.dim != other.dim:
            raise err_.MatrixDimensionsMismatch()

        return True if self * other.t == identity_matrix(self.n_rows, self.n_cols) else False

    def get_numpy_compatible_matrix(self, array_function):
        return array_function([i for i in self.elements])

    def is_positive_definite(self) -> bool:
        n_dimensions = self.n_rows

        if not self.is_symmetric():
            return False

        for i in range(1, n_dimensions + 1):
            minor = [row[:i] for row in self.elements[:i]]
            determinant_ = sum((minor[j][j] for j in range(i)))
            if determinant_ <= 0:
                return False

        return True

    def dot(self, other):
        if not isinstance(other, Matrix):
            raise ValueError('The other must be a Matrix object')

        return sum([i * j for i, j in zip(self.elements, other.elements)])


def determinant(matrix: Matrix or LList) -> float:
    def calculate_determinant(mat):
        n = len(mat)

        if n == 1:
            return mat[0][0]

        det = 0
        for j in range(n):
            sub_matrix = [row[:j] + row[j + 1:] for row in mat[1:]]
            cofactor = mat[0][j] * ((-1)**j)
            det += cofactor * calculate_determinant(sub_matrix)

        return det

    matrix = matrix.elements if isinstance(matrix, Matrix) else matrix

    n_rows = len(matrix)
    n_cols = len(matrix[0]) if matrix else 0

    if n_rows != n_cols:
        raise err_.NotASquareMatrix("Matrix must be square for determinant calculation.")

    return calculate_determinant(matrix)


class INVERSE:

    def __init__(self, matrix_elements: LList):
        self.elements = matrix_elements

    @staticmethod
    def _pivot_row(matrix, col):
        """Find the row with the largest absolute value in the current column."""
        max_val = 0
        max_row = -1
        for i in range(col, len(matrix)):
            if abs(matrix[i][col]) > max_val:
                max_val = abs(matrix[i][col])
                max_row = i
        return max_row

    @staticmethod
    def _swap_rows(matrix, i, j):
        """Swap two rows in the matrix."""
        matrix[i], matrix[j] = matrix[j], matrix[i]

    @staticmethod
    def _scale_row(matrix, row, factor):
        """Multiply a row by a scalar factor."""
        matrix[row] = [x * factor for x in matrix[row]]

    @staticmethod
    def _add_scaled_row(matrix, source_row, target_row, factor):
        """Add a scaled row to another row."""
        scaled_row = [x * factor for x in matrix[source_row]]
        matrix[target_row] = [x + y for x, y in zip(matrix[target_row], scaled_row)]

    def inverse(self):
        n_rows = len(self.elements)
        n_cols = len(self.elements[0])

        if n_rows != n_cols:
            raise err_.NotASquareMatrix("Matrix must be square for inverse calculation.")

        augmented_matrix = [row[:] + [int(i == j) for j in range(n_rows)] for i, row in enumerate(self.elements)]

        for col in range(n_rows):
            pivot_row = self._pivot_row(augmented_matrix, col)
            if pivot_row == -1:
                raise ValueError("Matrix is singular (no unique inverse).")

            self._swap_rows(augmented_matrix, col, pivot_row)

            pivot_value = augmented_matrix[col][col]
            self._scale_row(augmented_matrix, col, 1.0 / pivot_value)

            for row in range(n_rows):
                if row != col:
                    factor = -augmented_matrix[row][col]
                    self._add_scaled_row(augmented_matrix, col, row, factor)

        inverse_matrix = [row[n_rows:] for row in augmented_matrix]

        return inverse_matrix


class FRACTION:
    def __init__(self, decimal):
        self.fraction = Fraction(decimal).limit_denominator(100)

    def __str__(self):
        return str(self.fraction)

    @property
    def numerator(self):
        return self.fraction.numerator

    @property
    def denominator(self):
        return self.fraction.denominator


def identity_matrix(n_rows, n_cols):
    identity_ = Matrix([[0] * n_cols for _ in range(n_rows)])

    for i in range(n_rows):
        identity_[i][i] = 1

    return identity_


def null_matrix(n_rows, n_cols):
    mat_ = Matrix([[0] * n_cols for _ in range(n_rows)])
    if n_rows > 1:
        return mat_
    else:
        return mat_[0]


def vector_norm(vector: Matrix, squared: bool = False):
    if vector.n_cols == 1:
        import itertools
        vector_ = list(itertools.chain.from_iterable(vector.elements))
    else:
        vector_ = vector

    vec_norm_ = sum(i**2 for i in vector_.elements)

    return vec_norm_ if squared else sqrt(vec_norm_)
