# Interpolators

Interpolator methods are used to interpolate and find the value at a specified points given a set of $\langle x, y(x)
\rangle$.

| x   | A    | B    | C    | D    | E    | F    |
| --- | ---- | ---- | ---- | ---- | ---- | ---- |
| y   | y(A) | y(B) | y(C) | y(D) | y(E) | y(F) |

This library implements the following interpolators,

1. [Newton interpolators](https://github.com/syedalimohsinbukhari/Solvers/blob/md-docs/docs/interpolators/newton-interpolators.md)
2. [Lagrange interpolators](https://github.com/syedalimohsinbukhari/Solvers/blob/md-docs/docs/interpolators/lagrange-interpolators.md)
3. [Spline interpolators](https://github.com/syedalimohsinbukhari/Solvers/blob/md-docs/docs/interpolators/spline-interpolators.md)

## Using with x and y value tables

```python
from num_solvers.interpolators import newton_forward_interpolation

x_values = [0, 1, 2, 3, 4, 5]
y_values = [0, 1, 4, 9, 16, 25]

value_to_approximate = 2.2

result = newton_forward_interpolation(x_values, value_to_approximate, function_values=y_values)
print(result)
```

## Using with x values and function

```python
from num_solvers.interpolators import newton_backward_interpolation

x_values = [0, 1, 2, 3, 4, 5]
value_to_approximate = 4.3

result = newton_backward_interpolation(x_values, value_to_approximate, lambda x: x**2)
print(result)
```