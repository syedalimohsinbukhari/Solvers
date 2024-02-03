# Newton-Cotes integrators

The Newton-Cotes integrators (`NCIs`) are a set of formulae for numerical integration provided equally spaced intervals.

For a given single-variable function $f(\theta)$, the `NCIs` assume that the value of function $f$ defined at intervals $[a, b]$ is known at $n+1$ equally spaced points: $a \leq x_0 < x_1 < \cdots < x_n \leq b$. In this case, the `NCIs` are defined as;

$$\displaystyle \int_a^b f(\theta) \approx \sum_{i=0}^n \omega_i f(\theta_i)$$

For a closed `NCI` formulae, we have $\theta_i = a + ih$ and $h = (b-a)/n$. The number $h$ is called the `step-size` and $\omega_i$ are the weights.

This library implements, the following closed-form `NCIs`,

* [Trapezoidal rule](https://mathworld.wolfram.com/TrapezoidalRule.html)
* [Simpson's 1/3 rule](https://mathworld.wolfram.com/SimpsonsRule.html)
* [Simpson's 3/8 rule](https://mathworld.wolfram.com/Simpsons38Rule.html)
* [Boole's Rule](https://mathworld.wolfram.com/BoolesRule.html)
* [Weddle's Rule](https://mathworld.wolfram.com/WeddlesRule.html)[^1]

> [!TIP]:
> The `NCIs` are implemented in this library as composite rules. The user needs to specify the number of composites (multiples of general algorithms) to be used. Default is 1, which implements the default algorithms.

## Implementation

For a definite, single-varaible integral,

$$\int_a^bf(\theta)\,\text{d}\theta$$

The [SingleVariableNewtonCotesIntegrators](https://github.com/syedalimohsinbukhari/Solvers/blob/master/src/num_solvers/__backend/newton_cotes_.py#L85) class provides the functionality of integration using `NCIs` such as,

## Trapezoidal rule

$$\int_a^b f(\theta)\,\text{d}\theta \approx \frac{h}{2}\sum_{k=1}^N \Bigg\{f(x_{k-1}) + f(x_k)\Bigg\}$$

## Simpson's 1/3 rule

$$\int_a^b f(\theta)\,\text{d}\theta \approx \frac{h}{3}\sum_{k=1}^{N/2}\Bigg\{f(x_{2i-2}) + 4f(x_{2i-1}) + f(x_{2i})\Bigg\}$$

## Simpson's 3/8 rule

$$\int_a^b f(\theta)\,\text{d}\theta \approx \frac{3h}{8}\sum_{k=1}^{N/3}\Bigg\{f(x_{3i-3}) + 3f(x_{3i-2}) + 3f(x_{3i-1}) + f(x_{3i})\Bigg\}$$

## Boole's rule

$$\int_a^b f(\theta)\,\text{d}\theta \approx \frac{2h}{45} \sum_{k=1}^{N/4} \Bigg\{7f(x_{4k-4}) + 32f(x_{4k-3}) + 12 f(x_{4k-2}) + 32f(x_{4k-1}) + 7f(x_{4k})\Bigg\}$$

## Weddle's rule

$$\int_a^b f(\theta)\,\text{d}\theta \approx \frac{3h}{10} \sum_{k=1}^{N/6} \Bigg\{f(x_{6k-6}) + 5f(x_{6k-5}) + f(x_{6k-4}) + 6f(x_{6k-3}) + f(x_{6k-2}) + 4f(x_{6k-1}) + f(x_{6k})\Bigg\}$$


[^1]: Weddle's Rule has two implementations, one that uses `41/140` as a constant, and other which omits this factor in favor of simpler formulae with `3/10` as the constant.


