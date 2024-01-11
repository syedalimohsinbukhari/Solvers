"""Newton-Cotes integrators

This module implements several Newton-Cotes integrators, such as,

- trapezoidal_rule:
    The trapezoidal rule integrator.
- simpson_13:
    Simpson 1/3 rule integrator.
- simpson_38:
    Simpson 3/8 rule integrator.
- boole:
    Bool rule integrator.
- weddle_41h140:
    Weddle rule integrator.
- weddle_3h10:
    Modified Weddle rule integrator.

Created on Nov 29 14:48:53 2023
"""

__all__ = ['trapezoid_rule', 'simpson_13', 'simpson_38', 'boole', 'weddle_41h140', 'weddle_3h10']

from custom_inherit import doc_inherit

from .. import DOC_STYLE, Func, IFloat
from ..__backend.newton_cotes_ import solve_helper


def trapezoid_rule(function: Func, x_0: IFloat, x_n: IFloat, composites: int = 1, get_table: bool = False) -> IFloat:
    """Perform numerical integration using the Trapezoidal Rule.
    
    Parameters
    ----------
    function :
        The function to be integrated.
    x_0 :
        The lower limit of integration.
    x_n :
        The upper limit of integration.
    composites :
        Number of sub-intervals for composite method. Default is 1.
    get_table :
        Whether to print the table of values. Default is False.

    Returns
    -------
        The numerical result of integration.
    """

    return solve_helper(function, x_0, x_n, composites, get_table, 'trapezoid')


@doc_inherit(trapezoid_rule, style=DOC_STYLE)
def simpson_13(function: Func, x_0: IFloat, x_n: IFloat, composites: int = 1, get_table: bool = False) -> IFloat:
    """Perform numerical integration using Simpson's 1/3 Rule."""

    return solve_helper(function, x_0, x_n, composites, get_table, 'simpson_13')


@doc_inherit(trapezoid_rule, style=DOC_STYLE)
def simpson_38(function: Func, x_0: IFloat, x_n: IFloat, composites: int = 1, get_table: bool = False) -> IFloat:
    """Perform numerical integration using Simpson's 3/8 Rule."""

    return solve_helper(function, x_0, x_n, composites, get_table, 'simpson_38')


@doc_inherit(trapezoid_rule, style=DOC_STYLE)
def boole(function: Func, x_0: IFloat, x_n: IFloat, composites: int = 1, get_table: bool = False) -> IFloat:
    """Perform numerical integration using Boole's Rule."""

    return solve_helper(function, x_0, x_n, composites, get_table, 'boole')


@doc_inherit(trapezoid_rule, style=DOC_STYLE)
def weddle_41h140(function: Func, x_0: IFloat, x_n: IFloat, composites: int = 1, get_table: bool = False) -> IFloat:
    """Perform numerical integration using Weddle 41/140 Rule."""

    return solve_helper(function, x_0, x_n, composites, get_table, 'weddle_41h140')


@doc_inherit(trapezoid_rule, style=DOC_STYLE)
def weddle_3h10(function: Func, x_0: IFloat, x_n: IFloat, composites: int = 1, get_table: bool = False) -> IFloat:
    """Perform numerical integration using modified Weddle 3/10 Rule."""

    return solve_helper(function, x_0, x_n, composites, get_table, 'weddle_3h10')
