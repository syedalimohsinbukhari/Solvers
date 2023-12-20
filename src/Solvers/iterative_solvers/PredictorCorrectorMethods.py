"""Created on Dec 20 10:55:46 2023"""

import numpy as np

from .backend.PREDICTOR_CORRECTOR import get_x, trapezoidal_rule


def euler_method(ode, x_0, y_0, h=0.1, n=1, x_max=None, n_decimal: int = 6, full_result=False):
    h = (x_max - x_0) / n if h is None else h
    x_max = x_0 + n * h if x_max is None else x_max

    y_n, x_n = [y_0], [x_0]
    while x_0 < x_max:
        y_n.append(y_n[-1] + h * ode(x_n[-1], y_n[-1]))
        x_n.append(x_n[-1] + h)

        x_0 = np.round(x_0 + h, n_decimal)

    return (x_n, y_n) if full_result else y_n[-1]


def euler_trapezoidal(ode, x_0, y_0, h, n, n_decimal=6, get_predictor=False):
    x_n = get_x(x_0, n, h, n_decimal)
    y_n = np.zeros_like(x_n)

    predictor, corrector = np.zeros_like(x_n), np.zeros_like(x_n)
    y_n[0], predictor[0], corrector[0] = [y_0] * 3

    for i in range(1, n + 1):
        predictor_ = euler_method(ode=ode, x_0=x_n[i - 1], y_0=corrector[i - 1], h=h, n_decimal=n_decimal)
        predictor[i] = np.round(predictor_, n_decimal)

        corrector[i] = trapezoidal_rule(ode, x_n[i - 1], y_n[i - 1], x_n[i], predictor[i], h, n_decimal)
        y_n[i] = corrector[i]

    return (x_n, corrector) if not get_predictor else (x_n, (predictor, corrector))
