# [Runge-Kutta solvers](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods)

The Runge-Kutta (RK) solvers are a class of iterative solvers that are used to approximate solution of simultaneous non-linear equations.

The most common method from RK family is the Runge-Kutta method of fourth order (RK4). [Euler's method](https://en.wikipedia.org/wiki/Euler_method) is the most basic `explicit`, and simplest RK method.

This library has implementations of RK method for second, third and fourth order. These implementations are also extended to be applicable on any number of equations.

## RK-2

The RK2 method uses expansion of Taylor series up to order 2.

For an ODE, $y' = F(x, y)$, given that $y(x_0) = y_0$ the RK2 method can be written as,

$$ y_{n+1} = y_n + \frac{1}{2}\left(k_1 + k_2\right)$$

where,

$$ k_1 = hF\left(x_n, y_n\right)$$

$$ k_2 = hF\left(x_n + h, y_n + k_1\right)$$

## RK-3

Using Taylor expansion up to order 3, the RK-3 method can be written as,

$$ y_{n+1} = y_n + \frac{1}{6}\left(k1 + 4k_2 + k_3\right)$$

where,

$$ k_1 = hF\left(x_n, y_n\right)$$

$$ k_2 = hF\left(x_n + \frac{h}{2}, y_n + \frac{k_1}{2}\right)$$

$$ k_3 = hF\left(x_n + h, y_n - k_1 + 2k_2 \right)$$

## RK-4

Using Taylor expansion up to order 4, the RK-4 method can be written as,

$$ y_{n+1} = y_n + \frac{1}{6}\left(k1 + 2k_2 + 2k_3+k_4\right)$$

where,

$$ k_1 = hF\left(x_n, y_n\right)$$

$$ k_2 = hF\left(x_n + \frac{h}{2}, y_n + \frac{k_1}{2}\right)$$

$$ k_3 = hF\left(x_n + \frac{h}{2}, y_n + \frac{k_2}{2}\right)$$

$$ k_4 = hF\left(x_n + h, y_n + k_3 \right)$$

> [!NOTE]
> See also, [Multi-RK implementation](https://github.com/syedalimohsinbukhari/Solvers/blob/md-docs/docs/rk-multi-solvers.md)