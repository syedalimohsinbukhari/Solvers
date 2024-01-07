"""Lagrange's interpolation method

This module provides the functionality to implement the Lagrange Interpolation method via,

- lagrange_interpolation: A function to compute the interpolation value using Lagrange method.

Created on Nov 01 16:47:03 2023
"""

__all__ = ['lagrange_interpolation']

from .. import FList, IFloat, OptFunc, OptList
from ..__backend.interpolation_ import LagrangeInterpolation


def lagrange_interpolation(given_values: FList, value_to_approximate: IFloat, function: OptFunc = None,
                           function_values: OptList = None):

    interpolator = LagrangeInterpolation(given_values, value_to_approximate, function, function_values)

    return interpolator.interpolate()
