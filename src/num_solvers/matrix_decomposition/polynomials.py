"""Created on Jan 22 23:12:46 2024"""

from .matrix import InFractions
from .. import FList, IFloat, N_DECIMAL
from ..__backend.extra_ import round_list_, round_value_


class Polynomial:

    def __init__(self, polynomial: FList):
        self.polynomial = polynomial

    def __repr__(self):
        r_poly_ = round_list_(self.polynomial, 3)

        terms = [f"{coefficients:+}x^{deg_}"
                 if deg_ > 1 else f"{coefficients:+}x" if deg_ == 1 else f'{coefficients:+}'
                 for deg_, coefficients in enumerate(r_poly_[::-1])]

        polynomial_str = " ".join(terms[::-1])

        return polynomial_str

    @property
    def degree(self) -> IFloat:
        return len(self.polynomial) - 1

    @property
    def in_fractions(self):
        return [InFractions(j) for j in self.polynomial]

    @staticmethod
    def _give_output(output):
        return Polynomial(output)

    def _derivative(self):
        poly_deg = self.degree

        derivative_ = [(poly_deg - index) * coefficients for index, coefficients in enumerate(self.polynomial)]

        if derivative_[-1] == 0:
            del derivative_[-1]

        return self._give_output(derivative_)

    def _integrate(self):
        poly_deg = self.degree

        integral_ = [coefficients / (poly_deg - index + 1) for index, coefficients in enumerate(self.polynomial)]

        return self._give_output(integral_ + [0])

    def derivative(self):
        return self._derivative()

    def integrate(self):
        return self._integrate()

    def evaluate(self, value: IFloat, n_decimal: int = N_DECIMAL) -> IFloat:
        poly_degree = self.degree

        eval_ = [coefficient * value**(poly_degree - index) for index, coefficient in enumerate(self.polynomial[:-1])]

        return round_value_(sum(eval_ + [self.polynomial[-1]]), n_decimal)


c = Polynomial([17, 0, 566, 1, -2, -55, 4, -18])
d = c.integrate()
print(d.evaluate(2, 3))
print(d)
