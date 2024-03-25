"""Spline interpolation methods

This module provides the interpolation methods for linear, quadratic and natural cubic spline,
- linear_spline_interpolation
- quadratic_spline_interpolation
- natural_cubic_spline_interpolation

Created on Nov 01 23:03:42 2023
"""

__all__ = ['linear_spline_interpolation', 'quadratic_spline_interpolation', 'natural_cubic_spline_interpolation']

from custom_inherit import doc_inherit

from .. import DOC_STYLE, FList, IFloat, OptFunc, OptList
from ..__backend.spline_ import LinearSpline, NaturalCubicSpline, QuadraticSpline, SPLINE


@doc_inherit(SPLINE.__init__, style=DOC_STYLE)
def linear_spline_interpolation(given_values: FList, value_to_approximate: IFloat, function: OptFunc = None,
                                function_values: OptList = None, show_splines: bool = False) -> IFloat:
    """
    Method to perform linear spline interpolation.

    Parameters
    ----------

    show_splines: bool
        Whether to show the spline equations or not. Defaults to False.

    Returns
    -------
        The interpolated value.

    Examples
    --------
    1. Perform the linear spline interpolation for x^2 on x = [1, 2, 3, 4, 5], approximate 2.2.

    >>> linear_spline_interpolation(given_values=[1, 2, 3, 4, 5],
                                    value_to_approximate=2.2,
                                    function=lambda x: x**2)

    2. Perform the linear spline interpolation on x = [1, 2, 3, 4, 5], and y = [4, 6, 2, 4, 8, 10], approximate 3.3.

    >>> linear_spline_interpolation(given_values=[1, 2, 3, 4, 5],
                                    value_to_approximate=3.3,
                                    function_values=[4, 6, 2, 4, 8, 10])
    """

    linear_spline = LinearSpline(given_values, value_to_approximate, function, function_values)

    if show_splines:
        linear_spline.show_splines()

    return linear_spline.interpolate()


@doc_inherit(linear_spline_interpolation, style=DOC_STYLE)
def quadratic_spline_interpolation(given_values: FList, value_to_approximate: IFloat,
                                   function: OptFunc = None, function_values: OptList = None,
                                   last_equation: str = 'first', show_splines: bool = False) -> IFloat:
    """
    Perform quadratic spline interpolation.

    Examples
    --------
    1. Perform the quadratic spline interpolation for x^2 on x = [1, 2, 3, 4, 5], approximate 2.2.

    >>> quadratic_spline_interpolation(given_values=[1, 2, 3, 4, 5],
                                       value_to_approximate=2.2,
                                       function=lambda x: x**2)

    2. Perform the quadratic spline interpolation on x = [1, 2, 3, 4, 5], and y = [4, 6, 2, 4, 8, 10], approximate 3.3.

    >>> quadratic_spline_interpolation(given_values=[1, 2, 3, 4, 5],
                                       value_to_approximate=3.3,
                                       function_values=[4, 6, 2, 4, 8, 10])
    """

    quadratic_spline = QuadraticSpline(given_values, value_to_approximate, function, function_values, last_equation)

    if show_splines:
        quadratic_spline.show_splines()

    return quadratic_spline.interpolate()


@doc_inherit(linear_spline_interpolation, style=DOC_STYLE)
def natural_cubic_spline_interpolation(given_values: FList, value_to_approximate: IFloat,
                                       function: OptFunc = None, function_values: OptList = None,
                                       last_equation: str = 'first', show_splines: bool = False) -> IFloat:
    """
    Perform natural cubic spline interpolation.

    Examples
    --------
    1. Perform the natural cubic spline interpolation for x^2 on x = [1, 2, 3, 4, 5], approximate 2.2.

    >>> natural_cubic_spline_interpolation(given_values=[1, 2, 3, 4, 5],
                                           value_to_approximate=2.2,
                                           function=lambda x: x**2)

    2. Perform the natural cubic spline interpolation on x = [1, 2, 3, 4], and y = [4, 6, 2, 4, 8], approximate 3.3.

    >>> natural_cubic_spline_interpolation(given_values=[1, 2, 3, 4, 5],
                                           value_to_approximate=3.3,
                                           function_values=[4, 6, 2, 4, 8, 10])
    """

    natural_cubic_spline = NaturalCubicSpline(given_values, value_to_approximate, function, function_values,
                                              last_equation)

    if show_splines:
        natural_cubic_spline.show_splines()

    return natural_cubic_spline.interpolate()
