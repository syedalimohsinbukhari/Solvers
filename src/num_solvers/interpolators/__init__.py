"""Created on Dec 20 13:36:09 2023"""

from .lagrange_interpolation import lagrange_interpolation
from .newton_interpolation import newton_backward_interpolation, divided_difference_interpolation, newton_forward_interpolation
from .spline_interpolation import (linear_spline_interpolation, natural_cubic_spline_interpolation,
                                   quadratic_spline_interpolation)
