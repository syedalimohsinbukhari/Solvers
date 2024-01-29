"""Polynomials

This module provides functionality to create a polynomials and perform various operations on them. The base class is,

- :class:`Polynomial`

which gives the functionality of generating the polynomials. All polynomials have the following associated properties,

- degree: The degree of polynomial.
- in_fractions: Polynomial as a list of fractions.

Along with these properties, the polynomial object also has the following functions,

- derivative: Calculates the derivative of the polynomial.
- integral: Calculates the integral of the polynomial.
- is_zero: Checks whether the polynomial is 0 or not, e.g., all the coefficients of polynomial == 0.

The polynomial class can add two polynomial of different degrees,

Created on Jan 22 23:12:46 2024
"""

from fractions import Fraction

from .core_helpers_ import round_list_, round_value_
from .. import FList, IFloat, N_DECIMAL


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

    def __add__(self, other):
        if isinstance(other, Polynomial):
            if self.degree < other.degree:
                short_poly, long_poly = self, other
            else:
                short_poly, long_poly = other, self

            diff = long_poly.degree - short_poly.degree
            for _ in range(diff):
                short_poly.polynomial.insert(0, 0)

            return self._get_output([sum(x) for x in zip(self.polynomial, other.polynomial)])

        elif isinstance(other, (int, float)):
            temp_ = self.polynomial[:]
            temp_[-1] += other
            return temp_

    def __sub__(self, other):
        return self.__add__(-1 * other)

    def __truediv__(self, other):
        if isinstance(other, (int, float)):
            return self._get_output([coefficient / other for coefficient in self.polynomial])

    def __mul__(self, other):
        s_poly = self.polynomial

        if isinstance(other, Polynomial):
            o_poly = other.polynomial

            result_degree = len(s_poly) + len(o_poly) - 2
            result_coefficients = [0] * (result_degree + 1)

            for deg_self, c_self in enumerate(s_poly):
                for deg_other, c_other in enumerate(o_poly):
                    result_coefficients[deg_self + deg_other] += c_self * c_other

            return self._get_output(result_coefficients)

        if isinstance(other, (int, float)):
            return self._get_output([other * coefficient for coefficient in s_poly])

    __rmul__ = __mul__

    @staticmethod
    def _get_output(value):
        return Polynomial(value)

    @property
    def degree(self) -> IFloat:
        """
        Finds the degree of the given polynomial.
        Returns
        -------
            Degree of the polynomial.
        """

        return len(self.polynomial) - 1

    @property
    def in_fractions(self):
        """
        Displays the polynomial in fractions.
        Returns
        -------
            Polynomial in fractions.
        """

        return [InFractions(j) for j in self.polynomial]

    def _derivative(self):
        poly, poly_deg = self.polynomial, self.degree

        derivative_ = [(poly_deg - index) * coefficients for index, coefficients in enumerate(poly)]

        if derivative_[-1] == 0:
            del derivative_[-1]

        return self._get_output(derivative_)

    def _integral(self):
        poly, poly_deg = self.polynomial, self.degree

        integral_ = [coefficients / (poly_deg - index + 1) for index, coefficients in enumerate(poly)]

        return self._get_output(integral_ + [0])

    def is_zero(self):
        """Checks if the polynomial is zero polynomial."""
        return all(coefficients == 0 for coefficients in self.polynomial)

    def derivative(self):
        """Differentiates the given polynomial.

        Returns
        -------
            Differential of the given polynomial.
        """

        return self._derivative()

    def integral(self):
        """Integrates the given polynomial.

        Returns
        -------
            Integral of the given polynomial.
        """

        return self._integral()

    def evaluate(self, value: IFloat, n_decimal: int = N_DECIMAL) -> IFloat:
        """
        Evaluates the polynomial at the given value.

        Parameters
        ----------
        value:
            The value at which the polynomial should be evaluated.
        n_decimal:
            The number of decimal places to round off to.

        Returns
        -------
            Value of the polynomial after evaluation.
        """

        poly, poly_degree = self.polynomial, self.degree

        eval_ = [coefficient * value**(poly_degree - index) for index, coefficient in enumerate(poly[:-1])]

        return round_value_(sum(eval_ + [poly[-1]]), n_decimal)


class InFractions:
    def __init__(self, decimal_value: IFloat):
        self.fraction = Fraction(decimal_value).limit_denominator()

    def __repr__(self) -> str:
        return str(self.fraction)

    @property
    def numerator(self) -> IFloat:
        return self.fraction.numerator

    @property
    def denominator(self) -> IFloat:
        return self.fraction.denominator
