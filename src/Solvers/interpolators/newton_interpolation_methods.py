"""Created on Nov 02 00:25:14 2023"""

from custom_inherit import doc_inherit

from ._backend.interpolation_ import BkwInterpolation, DividedInterpolation, FwdInterpolation, Interpolation
from .. import DOC_STYLE, FList, IFloat, OptFunc, OptList


# TODO: Add module level docstrings

@doc_inherit(Interpolation.__init__, style=DOC_STYLE)
def forward_interpolation(given_values: FList, value_to_approximate: IFloat, function: OptFunc = None,
                          function_values: OptList = None, n_decimal: int = 8) -> IFloat:
    """
    Perform forward difference interpolation.

    Returns
    -------
        The interpolated result.
    """

    interpolation = FwdInterpolation(given_values, value_to_approximate, function, function_values, n_decimal)
    return interpolation.interpolate()


@doc_inherit(forward_interpolation, style=DOC_STYLE)
def backward_interpolation(given_values: FList, value_to_approximate: IFloat, function: OptFunc = None,
                           function_values: OptList = None):
    """Perform backwards difference interpolation."""

    interpolation = BkwInterpolation(given_values, value_to_approximate, function, function_values)
    return interpolation.interpolate()


@doc_inherit(forward_interpolation, style=DOC_STYLE)
def divided_interpolation(given_values: FList, value_to_approximate: IFloat, function: OptFunc = None,
                          function_values: OptList = None):
    """Perform divided difference interpolation."""

    interpolation = DividedInterpolation(given_values, value_to_approximate, function, function_values)
    return interpolation.interpolate()
