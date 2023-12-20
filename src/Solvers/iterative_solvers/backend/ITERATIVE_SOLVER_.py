"""Created on Oct 14 06:07:34 2023"""

from typing import List


class SysEqnSolver:
    def __init__(self, system_of_equations: List[List], solution: List, n_iter: int = 500, tol: float = 1e-5,
                 initial_guess: tuple = (0, 0, 0)):
        self.tol = tol
        self.n_iter = n_iter
        self.solution = solution
        self.initial_guess = initial_guess
        self.system_of_equations = system_of_equations
        self._arr1, self._arr2, self._arr3 = [self.initial_guess[0]], [self.initial_guess[1]], [self.initial_guess[2]]

    @property
    def solution_set(self):
        if not self._arr1 or len(self._arr1) == 1:
            self.solve()

        return self._arr1, self._arr2, self._arr3

    def _evaluate(self):
        equations, solution, initial_guess = self.system_of_equations, self.solution, self.initial_guess

        iter1_ = solution[0]
        iter1_ -= equations[0][1] * initial_guess[1]
        iter1_ -= equations[0][2] * initial_guess[2]
        iter1_ = iter1_ * equations[0][0]**-1

        iter2_ = solution[1]

        if self.__class__.__name__ == 'GaussJacobi':
            iter2_ -= equations[1][0] * initial_guess[0]
        else:
            iter2_ -= equations[1][0] * iter1_

        iter2_ -= equations[1][2] * initial_guess[2]
        iter2_ = iter2_ * equations[1][1]**-1

        iter3_ = solution[2]

        if self.__class__.__name__ == 'GaussJacobi':
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

    # def _string(self):
    #     equations, solution, initial_guess = self.system_of_equations, self.solution, self.initial_guess
    #     str1_ = f'({solution[0][0]}'
    #     str1_ += f' - ({equations[0][1]}*{initial_guess[1]})'
    #     str1_ += f' - ({equations[0][2]}*{initial_guess[2]}))'
    #
    #     str1_ = f'({Fraction(1, equations[0][0])})*{str1_}'
    #
    #     ic(str1_, eval(str1_))
    #
    #     str2_ = f'({solution[1][0]}'
    #
    #     if self.__class__.__name__ == 'GaussJacobi':
    #         str2_ += f' - ({equations[1][0]}*{initial_guess[0]})'
    #     else:
    #         str2_ += f' - ({equations[1][0]}*{eval(str1_)})'
    #
    #     str2_ += f' - ({equations[1][2]}*{initial_guess[2]}))'
    #
    #     str2_ = f'({Fraction(1, equations[1][1])})*{str2_}'
    #
    #     ic(str2_, eval(str2_))
    #
    #     str3_ = f'({solution[2][0]}'
    #
    #     if self.__class__.__name__ == 'GaussJacobi':
    #         str3_ += f' - ({equations[2][0]}*{initial_guess[0]})'
    #         str3_ += f' - ({equations[2][1]}*{initial_guess[1]}))'
    #     else:
    #         str3_ += f' - ({equations[2][0]}*{eval(str1_)})'
    #         str3_ += f' - ({equations[2][1]}*{eval(str2_)}))'
    #
    #     str3_ = f'({Fraction(1, equations[2][2])})*{str3_}'
    #
    #     ic(str3_, eval(str3_))

    def solve(self):
        for iter_ in range(self.n_iter):
            self.initial_guess = self._evaluate()

        # self.initial_guess = [round(i) for i in self.initial_guess]

        return self.initial_guess


class GaussJacobi(SysEqnSolver):
    pass


class GaussSeidel(SysEqnSolver):
    pass
