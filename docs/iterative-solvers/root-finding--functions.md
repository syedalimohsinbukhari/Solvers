# Root finding algorithms

In numerical analysis, a root-finding algorithm is an algorithm for finding zeros, also called "roots", of [continuous functions](https://en.wikipedia.org/wiki/Continuous_function). A [zero of a function](https://en.wikipedia.org/wiki/Zero_of_a_function) $f$, from the real numbers to real numbers or from the complex numbers to the complex numbers, is a number $x$ such that $f(x) = 0$.

Most numerical root-finding methods use iteration, producing a sequence of numbers that hopefully converges towards the root as its limit. They require one or more initial guesses of the root as starting values, then each iteration of the algorithm produces a successively more accurate approximation to the root. Since the iteration must be stopped at some point, these methods produce an approximation to the root, not an exact solution.

This library includes the following root-finding algorithms,

1. Bisection method
2. False position method (Regula-falsi method)
3. Secant method
4. Generalized secant method
5. Muller method
6. Newton-Raphson method
7. Ridder method
8. Sidi method
9. Steffensen method

## [Bisection method](https://en.wikipedia.org/wiki/Bisection_method)

1. Calculate $c$, the midpoint of the given interval, $c = \dfrac{a+b}{2}$.
2. Calculate the funciton value at midpoint, $f(c)$.
3. If $f(c) = 0$ or $\dfrac{b-a}{2} < \epsilon$, return $c$ as root.
4. If $f(c) < 0$, $a = c$ else $b = c$, repeat.

## [False-position method](https://en.wikipedia.org/wiki/Regula_falsi#The_regula_falsi_(false_position)_method)

The root of a given function $f(x)$ can be found using Regula-Falsi method as follows,

$$ c_k = \frac{a_kf(b_k) - b_kf(a_k)}{f(b_k) - f(a_k)} $$

The iterative procedure is as follows,
1. Initialize the intervals $[a, b]$ such that $f(a)\cdot f(b) < 0$
2. Repeat
   1. Estimate $c$ via the formula above,
   2. If $f(c) = 0$ or $\dfrac{b-a}{2} < \epsilon$, then $c$ is the root.
   3. Otherwise, if $f(a)\cdot f(c) < 0$, set $b = c$, else set $a = c$.

## [Seant method](https://en.wikipedia.org/wiki/Secant_method)

$$ x_{n+1} = x_n - f(x_n)\frac{x_n - x_{n-1}}{f(x_n) - f(x_{n-1})} $$

The iterative procedure is as follows,

Repeat
1. Estimate $x_{n+1}$ via the formula above,
2. Set, $x_{n-1} = x_n$ and $x_n = x_{n+1}$.
3. If $f(x_{n+1}) = 0$ or $f(x_{n+1}) < \epsilon$, $x_{n+1}$ is the root.

> [!NOTE]
> The secant method can also be written in terms of divided differnece
> $$x_{n+1} = x_n - \frac{f(x_n)}{f\left[x_n, x_{n-1}\right]}$$

## [Sidi's method](https://en.wikipedia.org/wiki/Sidi%27s_generalized_secant_method)

The Sidi's method is the generalization of the secand method, and can be iteratively performed as follows,

Repeat,
1. Calculate $x_{n+1}$ using the secand method,
2. Calculate $f(x_{n+1})$, and divided differences at $x_{n-1}$, $x_n$, $x_{n+1}$,
3. Calculate
$$x_{n+2} = x_{n+1} - \frac{f(x_{n+1})}{f[x_{n+1}, x_{n-1}] + (x_{n+1} - x_n)f[x_{n+1}, x_n, x_{n-1}]}$$
4. Set $x_{n-1} = x_n$, $x_n = x_{n+1}$, and $x_{n+1} = x_{n+2}$
5. if $f[x_{n+1}, x_{n-1}] + f[x_{n+1}, x_n, x_{n-1}] < \epsilon$, return $x_{n+2}$ as root.

## [Ridders' method](https://en.wikipedia.org/wiki/Ridders%27_method)

Repeat,
1. Calculate $x_\text{mid} = \dfrac{x_n + x_{n+1}}{2}$
2. Calculate,
$$ X = \dfrac{(x_\text{mid} - x_n) f(x_\text{mid})}{\sqrt{f(x_\text{mid})^2 - f(x_n)f(x_{n+1})}} $$
3. Evaluate
$$x_\text{new} = x_\text{mid} + \text{sign}\left[f(x_n), f(x_{n+1})\right]X$$
4. If, $f(x_\text{new}) = 0$ or $\sqrt{f(x_\text{mid})^2 - f(x_n)f(x_{n+1})} = 0$ return $x_\text{new}$ as the root.
5. Otherwise,
   1. If $f(x_\text{mid})f(x_\text{new}) < 0$, set $x_{n} = x_\text{mid}$ and $x_{n+1} = x_\text{new}$
   3. Else set $x_n = x_\text{new}$

> [!NOTE]
> $$ 
> \begin{equation*}
> \text{sign}\left[f(x_n), f(x_{n+1})\right] = 
>    \begin{cases}
>        + ,& f(x_n)f(x_{n+1}) > 0\\
>        - ,& \text{otherwise}
>    \end{cases}
> \end{equation*}
> $$

## Steffensen method

The Steffensen method is implemented in this library via two methods,

1. Root between two points
2. Root from initial guess

### Root between two points

$$ x_1 = x_0 - \frac{f(x_0)^2}{f\left(x_0 + f(x_0)\right) - f\left(x_0\right)} $$

1. Take two points $[a, b]$ such that, $a < b$ and $f(a)\cdot f(b) < 0$
2. Repeat,
   1. Calculate $x_0 = \dfrac{a+b}{2}$,
   2. Calculate $f(x_0)$ and $f(x_0 + f(x_0))$, and calculate $x_1$ using the above formula
   3. If $f(x_1) = 0$ than $x_1$ is an exact root, else $x_0 = x_1$.

### Root from initial guess