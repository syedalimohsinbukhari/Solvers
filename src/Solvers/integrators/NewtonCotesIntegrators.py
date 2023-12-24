"""Created on Nov 29 14:48:53 2023"""

from typing import Callable, List, Union

from ._backend.NEWTON_COTES import SingleVariableNewtonCotesIntegrators as SvNCI

INTEGRATOR_OUTPUT = Union[List[float], float]


def newton_raphson_solver(function: callable, derivative_of_function: callable, initial_guess: float,
                          tolerance: float = 1e-14) -> float:
    """
    Find the root of a function using the Newton-Raphson method.

    Parameters
    ----------
    function : callable
        The function for which the root is being sought.
    derivative_of_function : callable
        The derivative of the target function.
    initial_guess : float
        The initial guess for the root.
    tolerance : float, optional
        Tolerance for convergence. Default is 1e-14.

    Returns
    -------
    float
        The approximate root of the function.
    """
    f, df, x_0 = function, derivative_of_function, initial_guess

    if abs(f(x_0)) < tolerance:
        return x_0
    else:
        return newton_raphson_solver(f, df, x_0 - f(x_0) / df(x_0), tolerance)


def trapezoid_rule(function: Callable, x_0: float, x_n: float, composites: int = 1,
                   get_table: bool = False) -> INTEGRATOR_OUTPUT:
    """
    Perform numerical integration using the Trapezoidal Rule.

    Parameters
    ----------
    function : callable
        The function to be integrated.
    x_0 : float
        The lower limit of integration.
    x_n : float
        The upper limit of integration.
    composites : int, optional
        Number of sub-intervals for composite method. Default is 1.
    get_table : bool, optional
        Whether to return the table of values. Default is False.

    Returns
    -------
    INTEGRAL_OUTPUT
        The numerical result of integration or the table of values if ``get_table`` is True.
    """
    solver = SvNCI(function, x_0, x_n, composites, 'trapezoid')
    return solver.value_table() if get_table else solver.solve()


def simpson_13(function: Callable, x_0: float, x_n: float, composites: int = 1,
               get_table: bool = False) -> INTEGRATOR_OUTPUT:
    """
    Perform numerical integration using Simpson's 1/3 Rule.

    Parameters
    ----------
    function : callable
        The function to be integrated.
    x_0 : float
        The lower limit of integration.
    x_n : float
        The upper limit of integration.
    composites : int, optional
        Number of sub-intervals for composite method. Default is 1.
    get_table : bool, optional
        Whether to return the table of values. Default is False.

    Returns
    -------
    INTEGRAL_OUTPUT
        The numerical result of integration or the table of values if ``get_table`` is True.
    """
    solver = SvNCI(function, x_0, x_n, composites, 'simpson_13')
    return solver.value_table() if get_table else solver.solve()


def simpson_38(function: Callable, x_0: float, x_n: float, composites: int = 1,
               get_table: bool = False) -> INTEGRATOR_OUTPUT:
    """
    Perform numerical integration using Simpson's 3/8 Rule.

    Parameters
    ----------
    function : callable
        The function to be integrated.
    x_0 : float
        The lower limit of integration.
    x_n : float
        The upper limit of integration.
    composites : int, optional
        Number of sub-intervals for composite method. Default is 1.
    get_table : bool, optional
        Whether to return the table of values. Default is False.

    Returns
    -------
    INTEGRAL_OUTPUT
        The numerical result of integration or the table of values if ``get_table`` is True.
    """
    solver = SvNCI(function, x_0, x_n, composites, 'simpson_38')
    return solver.value_table() if get_table else solver.solve()


def boole(function: Callable, x_0: float, x_n: float, composites: int = 1,
          get_table: bool = False) -> INTEGRATOR_OUTPUT:
    """
    Perform numerical integration using Boole's Rule.

    Parameters
    ----------
    function : callable
        The function to be integrated.
    x_0 : float
        The lower limit of integration.
    x_n : float
        The upper limit of integration.
    composites : int, optional
        Number of sub-intervals for composite method. Default is 1.
    get_table : bool, optional
        Whether to return the table of values. Default is False.

    Returns
    -------
    INTEGRAL_OUTPUT
        The numerical result of integration or the table of values if ``get_table`` is True.
    """
    solver = SvNCI(function, x_0, x_n, composites, 'boole')
    return solver.value_table() if get_table else solver.solve()


def weddle_3h10(function: Callable, x_0: float, x_n: float, composites: int = 1,
                get_table: bool = False) -> INTEGRATOR_OUTPUT:
    """
    Perform numerical integration using Weddle's 3/10 Rule.

    Parameters
    ----------
    function : callable
        The function to be integrated.
    x_0 : float
        The lower limit of integration.
    x_n : float
        The upper limit of integration.
    composites : int, optional
        Number of sub-intervals for composite method. Default is 1.
    get_table : bool, optional
        Whether to return the table of values. Default is False.

    Returns
    -------
    INTEGRAL_OUTPUT
        The numerical result of integration or the table of values if ``get_table`` is True.
    """
    solver = SvNCI(function, x_0, x_n, composites, 'weddle_3h10')
    return solver.value_table() if get_table else solver.solve()


def weddle_41h140(function: Callable, x_0: float, x_n: float, composites: int = 1,
                  get_table: bool = False) -> INTEGRATOR_OUTPUT:
    """
    Perform numerical integration using Weddle's 41/140 Rule.

    Parameters
    ----------
    function : callable
        The function to be integrated.
    x_0 : float
        The lower limit of integration.
    x_n : float
        The upper limit of integration.
    composites : int, optional
        Number of sub-intervals for composite method. Default is 1.
    get_table : bool, optional
        Whether to return the table of values. Default is False.

    Returns
    -------
    INTEGRAL_OUTPUT
        The numerical result of integration or the table of values if ``get_table`` is True.
    """
    solver = SvNCI(function, x_0, x_n, composites, 'weddle_41h140')
    return solver.value_table() if get_table else solver.solve()
