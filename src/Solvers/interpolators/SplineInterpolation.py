"""Created on Nov 01 23:03:42 2023"""

from ._backend.SPLINE_ import LinearSpline, NaturalCubicSpline, QuadraticSpline


# TODO: Add docstring

def linear_spline_interpolation(given_values, value_to_approximate, function=None, function_values=None,
                                show_splines=False):
    linear_spline = LinearSpline(given_values, value_to_approximate, function, function_values)

    if show_splines:
        linear_spline.show_splines()

    return linear_spline.interpolate()


def quadratic_spline_interpolation(given_values, value_to_approximate, function=None, function_values=None,
                                   show_splines=False):
    quadratic_spline = QuadraticSpline(given_values, value_to_approximate, function, function_values)

    if show_splines:
        quadratic_spline.show_splines()

    return quadratic_spline.interpolate()


def natural_cubic_spline_interpolation(given_values, value_to_approximate, function=None, function_values=None,
                                       show_splines=False):
    natural_cubic_spline = NaturalCubicSpline(given_values, value_to_approximate, function, function_values)

    if show_splines:
        natural_cubic_spline.show_splines()

    return natural_cubic_spline.interpolate()
