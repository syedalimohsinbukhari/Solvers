# Predictor-Corrector method

In numerical analysis, predictor–corrector methods belong to a class of algorithms designed to integrate ordinary differential equations – to find an unknown function that satisfies a given differential equation. All such algorithms proceed in two steps:

1. The initial, *prediction* step, starts from a function fitted to the function-values and derivative-values at a preceding set of points to extrapolate ("anticipate") this function's value at a subsequent, new point.
2. The next, *corrector* step refines the initial approximation by using the predicted value of the function and another method to interpolate that unknown function's value at the same subsequent point.

This library implements an Euler-predictor and Trapezoidal-corrector method

## Euler-Trapezoidal PC method

Also known as [Heun's method](https://en.wikipedia.org/wiki/Heun%27s_method), this method uses the [Euler method](https://en.wikipedia.org/wiki/Euler_method) as a predictor, and [Trapezoidal rule](https://en.wikipedia.org/wiki/Trapezoidal_rule_(differential_equations)) as its corrector function.

Consider the differential equation,

$$ y^\prime = f(x, y),\qquad y(t_0) = y_0$$

The predictor step using Euler method can be obtained by

$$ Y^{\text{pred}}_{i+1} = y_i + hf(x_i, y_i)$$

Now we calculate the corrector step with Trapezoidal rule using the predictor step as,

$$
Y_{i+1}^{\text{corrector}} = y_i + \frac{h}{2}\left[f(x_i, y_i) + f(x_{i+1}, Y^{\text{pred}}_{i+1})\right]
$$

which is used as the first iteration $y_{i+1} = Y^{\text{corrector}}_{i+1}$.