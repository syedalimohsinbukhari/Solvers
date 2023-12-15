"""Created on Dec 16 03:28:44 2023"""


def rk2_solver(ode, x_0, y_0, h=0.1, x_max=1.0):
    y_n, x_n = [y_0], [x_0]

    while x_n[-1] + h < x_max:
        x_i, y_i = x_n[-1], y_n[-1]

        k1 = h * ode(x_i, y_i)
        k2 = h * ode(x_i + h, y_i + k1)

        y_n.append(y_i + 0.5 * (k1 + k2))
        x_n.append(x_i + h)

    return x_n, y_n


def rk3_solver(ode, x_0, y_0, h=0.1, x_max=1.0):
    y_n, x_n = [y_0], [x_0]

    while x_n[-1] + h < x_max:
        x_i, y_i = x_n[-1], y_n[-1]

        k1 = h * ode(x_i, y_i)
        k2 = h * ode(x_i + h / 2, y_i + k1 / 2)
        k3 = h * ode(x_i + h, y_i - k1 + 2 * k2)

        y_n.append(y_i + (1 / 6) * (k1 + 4 * k2 + k3))
        x_n.append(x_i + h)

    return x_n, y_n


def rk4_solver(ode, x_0, y_0, h=0.1, x_max=1.0):
    y_n, x_n = [y_0], [x_0]

    while x_n[-1] + h < x_max:
        x_i, y_i = x_n[-1], y_n[-1]

        k1 = h * ode(x_i, y_i)
        k2 = h * ode(x_i + h / 2, y_i + k1 / 2)
        k3 = h * ode(x_i + h / 2, y_i + k2 / 2)
        k4 = h * ode(x_i + h, y_i + k3)

        y_n.append(y_i + (1 / 6) * (k1 + 2 * k2 + 2 * k3 + k4))
        x_n.append(x_i + h)

    return x_n, y_n
