"""Lagrange's interpolation method

This module provides the functionality to implement the Lagrange Interpolation method via,

- lagrange_interpolation: A function to compute the interpolation value using Lagrange method.

Created on Nov 01 16:47:03 2023
"""

__all__ = ['lagrange_interpolation']

from custom_inherit import doc_inherit

from .. import DOC_STYLE, FList, IFloat, N_DECIMAL, OptFunc, OptList
from ..__backend.interpolation_ import Interpolation, LagrangeInterpolation


@doc_inherit(Interpolation.__init__, style=DOC_STYLE)
def lagrange_interpolation(given_values: FList, value_to_approximate: IFloat, function: OptFunc = None,
                           function_values: OptList = None, n_decimal: int = N_DECIMAL):
    """Performs interpolation via lagrange method.

    Returns
    -------
        The interpolated result.
    """

    interpolator = LagrangeInterpolation(given_values, value_to_approximate, function, function_values, n_decimal)

    return interpolator.interpolate()
