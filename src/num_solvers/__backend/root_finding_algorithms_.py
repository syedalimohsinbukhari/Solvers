"""Root finding algorithm - others

This module provides the divided difference backend method for root-finding algorithms.

- div_diff

Created on Jan 09 01:03:46 2024
"""

__all__ = ['div_diff']

from .. import FList, Func, IFloatOrFList


# Taken form https://en.wikipedia.org/wiki/Muller%27s_method#Computational_example
# minor tweaking applied on variable namings for consistency
def div_diff(function: Func, xs_: FList) -> IFloatOrFList:
    """Calculate the divided difference f[x0, x1, ...]."""
    if len(xs_) == 2:
        a, b = xs_
        return (function(a) - function(b)) / (a - b)
    else:
        f1 = div_diff(function, xs_[1:])
        f2 = div_diff(function, xs_[0:-1])

        return (f1 - f2) / (xs_[-1] - xs_[0])
