"""Created on Nov 29 14:48:53 2023"""

from .backend.NEWTON_COTES import NewtonCotesIterators


def newton_raphson_solver(function, func_derivative, x_0, tol=1e-14):
    f, df = function, func_derivative

    if abs(f(x_0)) < tol:
        return x_0
    else:
        return newton_raphson_solver(f, df, x_0 - f(x_0) / df(x_0), tol)


def trapezoid_rule(function, x_0, x_n, composites=1, get_table=False):
    solver = NewtonCotesIterators(function, x_0, x_n, composites, 'trapezoid')
    return solver.value_table() if get_table else solver.solve()


def simpson_13(function, x_0, x_n, composites=1, get_table=False):
    solver = NewtonCotesIterators(function, x_0, x_n, composites, 'simpson_13')
    return solver.value_table() if get_table else solver.solve()


def simpson_38(function, x_0, x_n, composites=1, get_table=False):
    solver = NewtonCotesIterators(function, x_0, x_n, composites, 'simpson_38')
    return solver.value_table() if get_table else solver.solve()


def boole(function, x_0, x_n, composites=1, get_table=False):
    solver = NewtonCotesIterators(function, x_0, x_n, composites, 'boole')
    return solver.value_table() if get_table else solver.solve()


def weddle_3h10(function, x_0, x_n, composites=1, get_table=False):
    solver = NewtonCotesIterators(function, x_0, x_n, composites, 'weddle_3h10')
    return solver.value_table() if get_table else solver.solve()


def weddle_41h140(function, x_0, x_n, composites=1, get_table=False):
    solver = NewtonCotesIterators(function, x_0, x_n, composites, 'weddle_41h140')
    return solver.value_table() if get_table else solver.solve()
