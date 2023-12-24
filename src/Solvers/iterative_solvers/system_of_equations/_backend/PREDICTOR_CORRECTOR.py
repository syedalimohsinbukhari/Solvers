"""Created on Dec 20 12:59:39 2023"""

import numpy as np


def _trapezoidal(ode, x_, y_, x_n1, y_n1, h):
    f_n = ode(x_, y_)
    f_p = ode(x_n1, y_n1)
    return y_ + (h / 2) * (f_n + f_p)


def trapezoidal_rule(ode, x_, y_, x_n1, y_n1, h, n_decimal):
    result = _trapezoidal(ode, x_, y_, x_n1, y_n1, h)
    return np.round(result, n_decimal)


def get_x(x_0, n, h, n_decimal):
    def np_round(x_):
        return np.round(x_, n_decimal)

    x_n = map(np_round, np.linspace(x_0, x_0 + n * h, n + 1))
    x_n = np.array(list(x_n))
    return np_round(x_n)
